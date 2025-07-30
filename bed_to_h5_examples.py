#!/usr/bin/env python3
"""
Example usage of the PLINK .bed to HDF5 converter.

This script demonstrates how to use the BedToH5Converter class
to convert PLINK files to HDF5 format programmatically.
"""

import os
import sys
from pathlib import Path

# Add the gensim directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from convert_bed_to_h5 import BedToH5Converter


def example_conversion():
    """Example of converting PLINK files to HDF5."""
    
    # Example 1: Convert files in the data directory (auto mode - recommended)
    input_dir = "data"
    output_dir = "data/h5_files"
    
    print("Example 1: Converting PLINK files to HDF5 (auto mode)")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    
    try:
        # Initialize converter
        converter = BedToH5Converter(
            input_dir=input_dir,
            output_dir=output_dir
        )
        
        # Find .bed files
        bed_files = converter.find_bed_files()
        print(f"\nFound {len(bed_files)} .bed files:")
        for bed_file, prefix in bed_files:
            print(f"  - {bed_file}")
        
        if bed_files:
            # Convert all files using auto mode
            print("\nStarting conversion (auto mode)...")
            converted_files = converter.convert_all_auto()
            
            print(f"\nConversion complete!")
            print(f"Converted {len(converted_files)} files:")
            for h5_file in converted_files:
                print(f"  - {h5_file}")
        else:
            print("No .bed files found for conversion")
    
    except Exception as e:
        print(f"Error during conversion: {e}")


def example_manual_modes():
    """Example of using specific processing modes."""
    
    print("\n" + "="*50)
    print("Example 1b: Manual processing mode selection")
    
    input_dir = "data"
    output_dir = "data/h5_files_manual"
    
    try:
        # Initialize converter
        converter = BedToH5Converter(
            input_dir=input_dir,
            output_dir=output_dir
        )
        
        # Find .bed files
        bed_files = converter.find_bed_files()
        
        if bed_files:
            print(f"Converting {len(bed_files)} files using standard processing...")
            converted_files = converter.convert_all()
            
            print(f"\nStandard conversion complete!")
            print(f"Converted {len(converted_files)} files")
        else:
            print("No .bed files found for manual conversion")
    
    except Exception as e:
        print(f"Error during manual conversion: {e}")


def example_chunked_conversion():
    """Example of chunked conversion for large datasets."""
    
    print("\n" + "="*50)
    print("Example 1b: Chunked conversion for large files")
    
    input_dir = "data"
    output_dir = "data/h5_files_chunked"
    
    try:
        # Initialize converter with chunked processing
        converter = BedToH5Converter(
            input_dir=input_dir,
            output_dir=output_dir,
            chunk_size=1000,  # Process 1000 samples at a time
            memory_usage=2.0   # Target 2GB memory usage
        )
        
        # Find .bed files
        bed_files = converter.find_bed_files()
        
        if bed_files:
            print(f"Converting {len(bed_files)} files using chunked processing...")
            converted_files = converter.convert_all_chunked()
            
            print(f"\nChunked conversion complete!")
            print(f"Converted {len(converted_files)} files:")
            for h5_file in converted_files:
                print(f"  - {h5_file}")
        else:
            print("No .bed files found for chunked conversion")
    
    except Exception as e:
        print(f"Error during chunked conversion: {e}")


def example_reading_h5():
    """Example of reading converted HDF5 files."""
    
    try:
        import h5py
        import numpy as np
    except ImportError:
        print("h5py is required for this example. Install with: pip install h5py")
        return
    
    print("\n" + "="*50)
    print("Example 2: Reading HDF5 files")
    
    # Find H5 files
    h5_files = list(Path("data").rglob("*.h5"))
    
    if not h5_files:
        print("No H5 files found. Run conversion first.")
        return
    
    # Read the first H5 file as an example
    h5_file = h5_files[0]
    print(f"\nReading H5 file: {h5_file}")
    
    try:
        with h5py.File(h5_file, 'r') as h5f:
            # Print file structure
            print("\nFile structure:")
            print(f"  Datasets: {list(h5f.keys())}")
            
            # Print metadata
            print("\nMetadata:")
            for attr_name, attr_value in h5f.attrs.items():
                print(f"  {attr_name}: {attr_value}")
            
            # Load genotype data
            genotypes = h5f['genotypes'][:]
            print(f"\nGenotype matrix shape: {genotypes.shape}")
            print(f"Data type: {genotypes.dtype}")
            
            # Show variant information
            if 'variants' in h5f:
                variants_group = h5f['variants']
                print(f"\nVariant information:")
                for dataset_name in variants_group.keys():
                    dataset = variants_group[dataset_name]
                    print(f"  {dataset_name}: {dataset.shape}")
            
            # Show sample information
            if 'samples' in h5f:
                samples_group = h5f['samples']
                print(f"\nSample information:")
                for dataset_name in samples_group.keys():
                    dataset = samples_group[dataset_name]
                    print(f"  {dataset_name}: {dataset.shape}")
    
    except Exception as e:
        print(f"Error reading H5 file: {e}")


def example_verification():
    """Example of verifying HDF5 files."""
    
    print("\n" + "="*50)
    print("Example 3: Verifying HDF5 files")
    
    # Find H5 files
    h5_files = list(Path("data").rglob("*.h5"))
    
    if not h5_files:
        print("No H5 files found for verification")
        return
    
    try:
        converter = BedToH5Converter(input_dir="data")
        
        print(f"Found {len(h5_files)} H5 files to verify:")
        
        valid_count = 0
        for h5_file in h5_files:
            print(f"\nVerifying: {h5_file}")
            is_valid = converter.verify_h5_file(str(h5_file))
            if is_valid:
                valid_count += 1
                print("  ✓ Valid")
            else:
                print("  ✗ Invalid")
        
        print(f"\nVerification summary: {valid_count}/{len(h5_files)} files are valid")
    
    except Exception as e:
        print(f"Error during verification: {e}")


if __name__ == "__main__":
    print("PLINK to HDF5 Converter Examples")
    print("="*50)
    
    # Check if data directory exists
    if not os.path.exists("data"):
        print("Data directory not found. Please ensure you have PLINK .bed files in a 'data' directory.")
        sys.exit(1)
    
    # Run examples
    example_conversion()
    example_manual_modes()
    example_chunked_conversion()
    example_reading_h5()
    example_verification()
    
    print("\n" + "="*50)
    print("Examples complete!")
    print("\nTo use the converter from command line:")
    print("  # Auto mode (recommended):")
    print("  python convert_bed_to_h5.py --input-dir data --output-dir data/h5_files --auto")
    print("  # Standard mode:")
    print("  python convert_bed_to_h5.py --input-dir data --output-dir data/h5_files")
    print("  # Chunked mode (for large files):")
    print("  python convert_bed_to_h5.py --input-dir data --output-dir data/h5_files --chunked --memory-usage 4.0")
    print("\nFor more options:")
    print("  python convert_bed_to_h5.py --help")
