#!/usr/bin/env python3
"""
Script to convert PLINK .bed files to HDF5 format.

This script recursively searches through a given directory for .bed files
(excluding temporary files) and converts them to HDF5 format for efficient
data storage and access.

Usage:
    python convert_bed_to_h5.py --input-dir /path/to/data --output-dir /path/to/output

Requirements:
    - pandas-plink
    - h5py
    - numpy
    - pandas
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Optional
import numpy as np
import pandas as pd

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

try:
    import h5py
except ImportError:
    print("Error: h5py is required. Install with: pip install h5py")
    sys.exit(1)

try:
    from pysnptools.snpreader import Bed
    PYSNPTOOLS_AVAILABLE = True
except ImportError:
    print("Error: pysnptools is required. Install with: pip install pysnptools")
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print("Warning: tqdm not available. Progress bars will be disabled.")
    tqdm = lambda x, **kwargs: x


class BedToH5Converter:
    """Convert PLINK .bed files to HDF5 format."""
    
    def __init__(self, input_dir: str, output_dir: Optional[str] = None, 
                 logger: Optional[logging.Logger] = None, chunk_size: Optional[int] = None,
                 memory_usage: float = 2.0):
        """
        Initialize the converter.
        
        Args:
            input_dir: Root directory to search for .bed files
            output_dir: Output directory for H5 files (default: same as input)
            logger: Optional logger instance
            chunk_size: Number of samples to process at once (auto-calculated if None)
            memory_usage: Target memory usage in GB for auto chunk size calculation
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir) if output_dir else self.input_dir
        self.logger = logger or self._setup_logger()
        self.chunk_size = chunk_size
        self.memory_usage = memory_usage
        
        if not self.input_dir.exists():
            raise ValueError(f"Input directory does not exist: {input_dir}")
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _setup_logger(self) -> logging.Logger:
        """Set up logging for the converter."""
        logger = logging.getLogger("bed_to_h5_converter")
        logger.setLevel(logging.INFO)
        
        # Clear existing handlers
        logger.handlers.clear()
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        return logger
    
    def calculate_chunk_size(self, n_samples: int, n_variants: int) -> int:
        """
        Calculate optimal chunk size based on memory usage.
        
        Args:
            n_samples: Number of samples in the dataset
            n_variants: Number of variants in the dataset
            
        Returns:
            Optimal chunk size for processing
        """
        if self.chunk_size is not None:
            return min(self.chunk_size, n_samples)
        
        # Estimate memory per sample (4 bytes per genotype * number of variants)
        bytes_per_sample = n_variants * 4  # float32
        
        # Target memory usage in bytes
        target_memory_bytes = self.memory_usage * 1024**3  # Convert GB to bytes
        
        # Calculate chunk size
        chunk_size = max(1, int(target_memory_bytes / bytes_per_sample))
        
        # Don't exceed total number of samples
        chunk_size = min(chunk_size, n_samples)
        
        # For very large datasets, be even more conservative
        if n_variants > 500000:  # More than 500K variants
            # Use maximum 10GB chunks for very large datasets
            max_chunk_memory_bytes = 10 * 1024**3  # 10GB
            max_chunk_size = max(1, int(max_chunk_memory_bytes / bytes_per_sample))
            chunk_size = min(chunk_size, max_chunk_size)
            self.logger.info(f"Large dataset detected ({n_variants:,} variants) - using conservative chunk size")
        
        self.logger.info(f"Calculated chunk size: {chunk_size:,} samples")
        self.logger.info(f"Memory per chunk: {(chunk_size * bytes_per_sample) / 1024**3:.2f} GB")
        
        return chunk_size
    
    def check_memory_requirements(self, n_samples: int, n_variants: int) -> Tuple[bool, float, float]:
        """
        Check if chunking is needed based on available memory and dataset size.
        
        Args:
            n_samples: Number of samples in the dataset
            n_variants: Number of variants in the dataset
            
        Returns:
            Tuple of (needs_chunking, required_memory_gb, available_memory_gb)
        """
        # Calculate memory required to load full dataset
        # Each genotype is float32 (4 bytes)
        required_memory_bytes = n_samples * n_variants * 4
        required_memory_gb = required_memory_bytes / (1024**3)
        
        # Check if we're running under SLURM and get allocated memory
        slurm_memory_gb = self._get_slurm_memory_limit()
        
        if slurm_memory_gb is not None:
            # We're running under SLURM - use allocated memory
            available_memory_gb = slurm_memory_gb
            self.logger.info(f"SLURM job detected: allocated {slurm_memory_gb:.1f} GB memory")
            self.logger.info(f"Dataset memory requirement: {required_memory_gb:.1f} GB")
            
            # Use safety factor - don't use more than 70% of allocated memory
            safety_factor = 0.7
            safe_memory_gb = available_memory_gb * safety_factor
            
        elif PSUTIL_AVAILABLE:
            # Fall back to system memory detection
            try:
                memory_info = psutil.virtual_memory()
                available_memory_gb = memory_info.available / (1024**3)
                total_memory_gb = memory_info.total / (1024**3)
                
                self.logger.info(f"System memory: {total_memory_gb:.1f} GB total, {available_memory_gb:.1f} GB available")
                self.logger.info(f"Dataset memory requirement: {required_memory_gb:.1f} GB")
                
                # Use safety factor - don't use more than 70% of available memory
                safety_factor = 0.7
                safe_memory_gb = available_memory_gb * safety_factor
                
            except Exception as e:
                self.logger.warning(f"Could not determine system memory: {e}")
                # Default to conservative chunking
                needs_chunking = required_memory_gb > 1.0
                return needs_chunking, required_memory_gb, 0.0
        else:
            self.logger.warning("Neither SLURM nor psutil available - cannot determine memory limits")
            # If we can't determine memory, be conservative and use chunking for datasets >1GB
            needs_chunking = required_memory_gb > 1.0
            return needs_chunking, required_memory_gb, 0.0
        
        needs_chunking = required_memory_gb > safe_memory_gb
        
        if needs_chunking:
            self.logger.info(f"Chunking required: dataset needs {required_memory_gb:.1f} GB, "
                           f"but only {safe_memory_gb:.1f} GB safely available")
        else:
            self.logger.info(f"Chunking not required: dataset fits in available memory "
                           f"({required_memory_gb:.1f} GB needed, {safe_memory_gb:.1f} GB available)")
        
        return needs_chunking, required_memory_gb, available_memory_gb
    
    def _get_slurm_memory_limit(self) -> Optional[float]:
        """
        Get the memory limit from SLURM environment variables.
        
        Returns:
            Memory limit in GB if running under SLURM, None otherwise
        """
        # Check for direct memory allocation first
        if 'SLURM_MEM_PER_NODE' in os.environ:
            try:
                memory_mb = int(os.environ['SLURM_MEM_PER_NODE'])
                return memory_mb / 1024.0  # Convert to GB
            except (ValueError, KeyError) as e:
                self.logger.warning(f"Could not parse SLURM_MEM_PER_NODE: {e}")
        
        # Check for memory per CPU allocation
        if 'SLURM_MEM_PER_CPU' in os.environ and 'SLURM_CPUS_PER_TASK' in os.environ:
            try:
                mem_per_cpu_mb = int(os.environ['SLURM_MEM_PER_CPU'])
                cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
                total_memory_mb = mem_per_cpu_mb * cpus
                self.logger.info(f"SLURM memory calculation: {mem_per_cpu_mb}MB/CPU × {cpus} CPUs = {total_memory_mb}MB total")
                return total_memory_mb / 1024.0  # Convert to GB
            except (ValueError, KeyError) as e:
                self.logger.warning(f"Could not parse SLURM memory per CPU variables: {e}")
        
        # Check if we're in a SLURM job but couldn't determine memory
        if 'SLURM_JOB_ID' in os.environ:
            self.logger.warning("Running in SLURM but could not determine memory allocation")
            self.logger.warning("Available SLURM variables:")
            for key, value in os.environ.items():
                if key.startswith('SLURM_'):
                    self.logger.warning(f"  {key}={value}")
        
        return None
    
    def auto_determine_processing_mode(self, bed_file_path: Path) -> str:
        """
        Automatically determine whether to use chunked or standard processing.
        
        Args:
            bed_file_path: Path to the .bed file
            
        Returns:
            Processing mode: 'chunked' or 'standard'
        """
        try:
            # Read metadata to get dataset dimensions using pysnptools
            plink_prefix = str(bed_file_path.with_suffix(''))
            bed_reader = Bed(plink_prefix, count_A1=False)  # Don't load data yet, just metadata
            
            n_samples = bed_reader.iid_count
            n_variants = bed_reader.sid_count
            
            # Check memory requirements
            needs_chunking, required_gb, available_gb = self.check_memory_requirements(n_samples, n_variants)
            
            if needs_chunking:
                # Auto-adjust memory usage target based on available memory
                if available_gb > 0:
                    # Use 50% of available memory for chunking
                    self.memory_usage = min(self.memory_usage, available_gb * 0.5)
                    self.logger.info(f"Adjusted memory usage target to {self.memory_usage:.1f} GB")
                
                return 'chunked'
            else:
                return 'standard'
                
        except Exception as e:
            self.logger.warning(f"Could not determine processing mode automatically: {e}")
            # Default to standard if we can't determine
            return 'standard'
            return 'standard'
    
    def find_bed_files(self) -> List[Tuple[Path, str]]:
        """
        Find all .bed files in the input directory and subdirectories.
        Excludes files containing 'temporary' in their name.
        
        Returns:
            List of tuples (bed_file_path, prefix_without_extension)
        """
        bed_files = []
        
        for bed_file in self.input_dir.rglob("*.bed"):
            # Skip temporary files
            if "temporary" in bed_file.name.lower():
                self.logger.info(f"Skipping temporary file: {bed_file}")
                continue
            
            # Get the prefix (filename without .bed extension)
            prefix = bed_file.stem
            
            # Check if corresponding .bim and .fam files exist
            bim_file = bed_file.with_suffix('.bim')
            fam_file = bed_file.with_suffix('.fam')
            
            if not bim_file.exists():
                self.logger.warning(f"Missing .bim file for {bed_file}, skipping")
                continue
            
            if not fam_file.exists():
                self.logger.warning(f"Missing .fam file for {bed_file}, skipping")
                continue
            
            bed_files.append((bed_file, prefix))
            self.logger.info(f"Found complete PLINK fileset: {bed_file.parent / prefix}")
        
        return bed_files
    
    def convert_bed_to_h5(self, bed_file_path: Path, prefix: str) -> str:
        """
        Convert a single PLINK .bed file to HDF5 format.
        
        Args:
            bed_file_path: Path to the .bed file
            prefix: Prefix for the output file
            
        Returns:
            Path to the created H5 file
        """
        try:
            self.logger.info(f"Converting {bed_file_path} to HDF5...")
            
            # Read PLINK files using pysnptools
            plink_prefix = str(bed_file_path.with_suffix(''))
            bed_reader = Bed(plink_prefix)
            
            # Read all data
            snp_data = bed_reader.read()
            genotype_matrix = snp_data.val  # Get the genotype matrix
            
            # Get sample and variant information
            sample_info = snp_data.iid  # Individual IDs
            variant_info = snp_data.sid  # SNP IDs
            
            # Create output H5 file path
            # Maintain the relative directory structure
            rel_path = bed_file_path.parent.relative_to(self.input_dir)
            output_subdir = self.output_dir / rel_path
            output_subdir.mkdir(parents=True, exist_ok=True)
            
            h5_file_path = output_subdir / f"{prefix}.h5"
            
            # Save to HDF5
            with h5py.File(h5_file_path, 'w') as h5f:
                # Save genotype data
                h5f.create_dataset('genotype_data', data=genotype_matrix,
                                    dtype='float32',
                                    compression='gzip')
            
            self.logger.info(f"Successfully converted to {h5_file_path}")
            self.logger.info(f"  - Samples: {genotype_matrix.shape[0]}")
            self.logger.info(f"  - Variants: {genotype_matrix.shape[1]}")
            
            return str(h5_file_path)
            
        except Exception as e:
            self.logger.error(f"Error converting {bed_file_path}: {e}")
            raise
    
    def convert_bed_to_h5_chunked(self, bed_file_path: Path, prefix: str) -> str:
        """
        Convert a single PLINK .bed file to HDF5 format using chunked processing.
        This approach is more memory-efficient for large datasets.
        
        Args:
            bed_file_path: Path to the .bed file
            prefix: Prefix for the output file
            
        Returns:
            Path to the created H5 file
        """
        try:
            self.logger.info(f"Converting {bed_file_path} to HDF5 (chunked)...")
            
            # Read PLINK file metadata using pysnptools
            plink_prefix = str(bed_file_path.with_suffix(''))
            bed_reader = Bed(plink_prefix)
            
            n_samples = bed_reader.iid_count
            n_variants = bed_reader.sid_count
            
            self.logger.info(f"Dataset dimensions: {n_samples:,} samples × {n_variants:,} variants")
            
            # Calculate optimal chunk size
            chunk_size = self.calculate_chunk_size(n_samples, n_variants)
            
            # Create output H5 file path
            rel_path = bed_file_path.parent.relative_to(self.input_dir)
            output_subdir = self.output_dir / rel_path
            output_subdir.mkdir(parents=True, exist_ok=True)
            
            h5_file_path = output_subdir / f"{prefix}.h5"
            
            # Process data in chunks and write to HDF5
            with h5py.File(h5_file_path, 'w') as h5f:
                # Initialize the dataset
                self.logger.info(f"Initializing HDF5 dataset with shape ({n_samples}, {n_variants})")
                dset = h5f.create_dataset('genotype_data',
                                        shape=(n_samples, n_variants),
                                        dtype='float32',
                                        compression='gzip',
                                        chunks=True)
                
                # Process samples in chunks
                for i in tqdm(range(0, n_samples, chunk_size), desc="Processing chunks"):
                    self.logger.info(f'Processing chunk {i//chunk_size + 1}: samples {i} to {min(i + chunk_size, n_samples)}')
                    
                    # Calculate end index for this chunk
                    end_idx = min(i + chunk_size, n_samples)
                    
                    # Read chunk of data using pysnptools slicing
                    chunk_data = bed_reader[i:end_idx, :].read().val
                    
                    # Write chunk data to the dataset
                    dset[i:end_idx, :] = chunk_data
            
            self.logger.info(f"Successfully converted to {h5_file_path}")
            self.logger.info(f"  - Samples: {n_samples:,}")
            self.logger.info(f"  - Variants: {n_variants:,}")
            self.logger.info(f"  - Processed in {(n_samples + chunk_size - 1) // chunk_size} chunks")
            
            return str(h5_file_path)
            
        except Exception as e:
            self.logger.error(f"Error converting {bed_file_path}: {e}")
            raise
    
    def convert_bed_to_h5_auto(self, bed_file_path: Path, prefix: str) -> str:
        """
        Convert a PLINK .bed file to HDF5 format using automatically determined processing mode.
        
        Args:
            bed_file_path: Path to the .bed file
            prefix: Prefix for the output file
            
        Returns:
            Path to the created H5 file
        """
        # Determine optimal processing mode
        processing_mode = self.auto_determine_processing_mode(bed_file_path)
        
        self.logger.info(f"Auto-selected processing mode: {processing_mode}")
        
        if processing_mode == 'chunked':
            return self.convert_bed_to_h5_chunked(bed_file_path, prefix)
        else:
            return self.convert_bed_to_h5(bed_file_path, prefix)
    
    def convert_all(self) -> List[str]:
        """
        Convert all found .bed files to HDF5 format.
        
        Returns:
            List of paths to created H5 files
        """
        bed_files = self.find_bed_files()
        
        if not bed_files:
            self.logger.warning("No .bed files found for conversion")
            return []
        
        self.logger.info(f"Found {len(bed_files)} .bed files to convert")
        
        converted_files = []
        failed_conversions = []
        
        for i, (bed_file_path, prefix) in enumerate(bed_files, 1):
            try:
                self.logger.info(f"Processing {i}/{len(bed_files)}: {bed_file_path}")
                h5_file = self.convert_bed_to_h5(bed_file_path, prefix)
                converted_files.append(h5_file)
                
            except Exception as e:
                self.logger.error(f"Failed to convert {bed_file_path}: {e}")
                failed_conversions.append(str(bed_file_path))
                continue
        
        # Summary
        self.logger.info(f"\nConversion complete!")
        self.logger.info(f"  - Successfully converted: {len(converted_files)}")
        self.logger.info(f"  - Failed conversions: {len(failed_conversions)}")
        
        if failed_conversions:
            self.logger.warning("Failed conversions:")
            for failed_file in failed_conversions:
                self.logger.warning(f"  - {failed_file}")
        
        return converted_files
    
    def convert_all_chunked(self) -> List[str]:
        """
        Convert all found .bed files to HDF5 format using chunked processing.
        
        Returns:
            List of paths to created H5 files
        """
        bed_files = self.find_bed_files()
        
        if not bed_files:
            self.logger.warning("No .bed files found for conversion")
            return []
        
        self.logger.info(f"Found {len(bed_files)} .bed files to convert (chunked mode)")
        
        converted_files = []
        failed_conversions = []
        
        for i, (bed_file_path, prefix) in enumerate(bed_files, 1):
            try:
                self.logger.info(f"Processing {i}/{len(bed_files)}: {bed_file_path}")
                h5_file = self.convert_bed_to_h5_chunked(bed_file_path, prefix)
                converted_files.append(h5_file)
                
            except Exception as e:
                self.logger.error(f"Failed to convert {bed_file_path}: {e}")
                failed_conversions.append(str(bed_file_path))
                continue
        
        # Summary
        self.logger.info(f"\nConversion complete!")
        self.logger.info(f"  - Successfully converted: {len(converted_files)}")
        self.logger.info(f"  - Failed conversions: {len(failed_conversions)}")
        
        if failed_conversions:
            self.logger.warning("Failed conversions:")
            for failed_file in failed_conversions:
                self.logger.warning(f"  - {failed_file}")
        
        return converted_files
    
    def convert_all_auto(self) -> List[str]:
        """
        Convert all found .bed files to HDF5 format using automatic processing mode selection.
        
        Returns:
            List of paths to created H5 files
        """
        bed_files = self.find_bed_files()
        
        if not bed_files:
            self.logger.warning("No .bed files found for conversion")
            return []
        
        self.logger.info(f"Found {len(bed_files)} .bed files to convert (auto mode)")
        
        converted_files = []
        failed_conversions = []
        
        for i, (bed_file_path, prefix) in enumerate(bed_files, 1):
            try:
                self.logger.info(f"Processing {i}/{len(bed_files)}: {bed_file_path}")
                h5_file = self.convert_bed_to_h5_auto(bed_file_path, prefix)
                converted_files.append(h5_file)
                
            except Exception as e:
                self.logger.error(f"Failed to convert {bed_file_path}: {e}")
                failed_conversions.append(str(bed_file_path))
                continue
        
        # Summary
        self.logger.info(f"\nConversion complete!")
        self.logger.info(f"  - Successfully converted: {len(converted_files)}")
        self.logger.info(f"  - Failed conversions: {len(failed_conversions)}")
        
        if failed_conversions:
            self.logger.warning("Failed conversions:")
            for failed_file in failed_conversions:
                self.logger.warning(f"  - {failed_file}")
        
        return converted_files
    
    def verify_h5_file(self, h5_file_path: str) -> bool:
        """
        Verify that an H5 file was created correctly.
        
        Args:
            h5_file_path: Path to the H5 file
            
        Returns:
            True if file is valid, False otherwise
        """
        try:
            with h5py.File(h5_file_path, 'r') as h5f:
                # Check required datasets exist
                required_datasets = ['genotype_data']
                
                for dataset in required_datasets:
                    if dataset not in h5f:
                        self.logger.error(f"Missing dataset '{dataset}' in {h5_file_path}")
                        return False

                self.logger.info(f"H5 file verification passed: {h5_file_path}")
                return True
                
        except Exception as e:
            self.logger.error(f"Error verifying {h5_file_path}: {e}")
            return False


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Convert PLINK .bed files to HDF5 format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert all .bed files in current directory
  python convert_bed_to_h5.py --input-dir .

  # Convert files and save to different directory
  python convert_bed_to_h5.py --input-dir data/simulations --output-dir data/h5_files

  # Auto mode - automatically determine processing method (recommended)
  python convert_bed_to_h5.py --input-dir data --auto

  # Use chunked processing for large files
  python convert_bed_to_h5.py --input-dir data --chunked --memory-usage 4.0

  # Use chunked processing with specific chunk size
  python convert_bed_to_h5.py --input-dir data --chunked --chunk-size 1000

  # Enable verbose logging
  python convert_bed_to_h5.py --input-dir data --verbose

  # Verify converted files
  python convert_bed_to_h5.py --input-dir data --verify-only
        """
    )
    
    parser.add_argument(
        '--input-dir', '-i',
        required=True,
        help='Input directory to search for .bed files'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        help='Output directory for H5 files (default: same as input)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--verify-only',
        action='store_true',
        help='Only verify existing H5 files, do not convert'
    )
    
    parser.add_argument(
        '--log-file',
        help='Log file path (default: console only)'
    )
    
    parser.add_argument(
        '--chunked',
        action='store_true',
        help='Use chunked processing for memory efficiency (recommended for large files)'
    )
    
    parser.add_argument(
        '--auto',
        action='store_true',
        help='Automatically determine processing mode based on available memory (recommended)'
    )
    
    parser.add_argument(
        '--chunk-size',
        type=int,
        help='Number of samples to process at once (auto-calculated if not specified)'
    )
    
    parser.add_argument(
        '--memory-usage',
        type=float,
        default=2.0,
        help='Target memory usage in GB for auto chunk size calculation (default: 2.0)'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = logging.getLogger("bed_to_h5_converter")
    level = logging.DEBUG if args.verbose else logging.INFO
    logger.setLevel(level)
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    try:
        # Initialize converter
        converter = BedToH5Converter(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            logger=logger,
            chunk_size=args.chunk_size,
            memory_usage=args.memory_usage
        )
        
        if args.verify_only:
            # Find and verify existing H5 files
            h5_files = list(converter.output_dir.rglob("*.h5"))
            logger.info(f"Found {len(h5_files)} H5 files to verify")
            
            valid_files = 0
            for h5_file in h5_files:
                if converter.verify_h5_file(str(h5_file)):
                    valid_files += 1
            
            logger.info(f"Verification complete: {valid_files}/{len(h5_files)} files valid")
            
        else:
            # Convert files
            if args.auto:
                logger.info("Using automatic processing mode (recommended)")
                converted_files = converter.convert_all_auto()
            elif args.chunked:
                logger.info("Using chunked processing mode")
                converted_files = converter.convert_all_chunked()
            else:
                logger.info("Using standard processing mode")
                converted_files = converter.convert_all()
            
            if converted_files:
                logger.info(f"\nConversion summary:")
                logger.info(f"Created {len(converted_files)} H5 files")
                
                # Optionally verify converted files
                logger.info("Verifying converted files...")
                verified = 0
                for h5_file in converted_files:
                    if converter.verify_h5_file(h5_file):
                        verified += 1
                
                logger.info(f"Verification: {verified}/{len(converted_files)} files valid")
    
    except KeyboardInterrupt:
        logger.info("Conversion interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
