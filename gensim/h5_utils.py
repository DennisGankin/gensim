"""
HDF5 utilities for reading converted PLINK files.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Any
from pathlib import Path

try:
    import h5py
    HDF5_AVAILABLE = True
except ImportError:
    HDF5_AVAILABLE = False


class H5PLINKReader:
    """Reader for HDF5 files created from PLINK data."""
    
    def __init__(self, h5_file_path: str):
        """
        Initialize the HDF5 PLINK reader.
        
        Args:
            h5_file_path: Path to the HDF5 file
        """
        if not HDF5_AVAILABLE:
            raise ImportError("h5py is required to read HDF5 files. Install with: pip install h5py")
        
        self.h5_file_path = Path(h5_file_path)
        if not self.h5_file_path.exists():
            raise FileNotFoundError(f"HDF5 file not found: {h5_file_path}")
        
        self._h5_file = None
        self._metadata = None
    
    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
    
    def open(self):
        """Open the HDF5 file."""
        if self._h5_file is None:
            self._h5_file = h5py.File(self.h5_file_path, 'r')
            self._load_metadata()
    
    def close(self):
        """Close the HDF5 file."""
        if self._h5_file is not None:
            self._h5_file.close()
            self._h5_file = None
    
    def _load_metadata(self):
        """Load metadata from the HDF5 file."""
        if self._h5_file is None:
            raise RuntimeError("HDF5 file is not open")
        
        self._metadata = dict(self._h5_file.attrs)
    
    @property
    def metadata(self) -> Dict[str, Any]:
        """Get file metadata."""
        if self._metadata is None:
            self._load_metadata()
        return self._metadata.copy()
    
    @property
    def n_samples(self) -> int:
        """Get number of samples."""
        return self.metadata.get('n_samples', 0)
    
    @property
    def n_variants(self) -> int:
        """Get number of variants."""
        return self.metadata.get('n_variants', 0)
    
    @property
    def shape(self) -> Tuple[int, int]:
        """Get genotype matrix shape (n_samples, n_variants)."""
        return (self.n_samples, self.n_variants)
    
    def get_genotypes(self, sample_indices: Optional[np.ndarray] = None,
                     variant_indices: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Get genotype data with optional subsetting.
        
        Args:
            sample_indices: Indices of samples to retrieve (all if None)
            variant_indices: Indices of variants to retrieve (all if None)
            
        Returns:
            Genotype matrix (n_samples, n_variants)
        """
        if self._h5_file is None:
            raise RuntimeError("HDF5 file is not open")
        
        genotypes = self._h5_file['genotypes']
        
        if sample_indices is None and variant_indices is None:
            return genotypes[:]
        elif sample_indices is None:
            return genotypes[:, variant_indices]
        elif variant_indices is None:
            return genotypes[sample_indices, :]
        else:
            return genotypes[np.ix_(sample_indices, variant_indices)]
    
    def get_variants_df(self) -> pd.DataFrame:
        """
        Get variant information as a pandas DataFrame.
        
        Returns:
            DataFrame with variant information (columns from .bim file)
        """
        if self._h5_file is None:
            raise RuntimeError("HDF5 file is not open")
        
        if 'variants' not in self._h5_file:
            raise ValueError("No variant information found in HDF5 file")
        
        variants_group = self._h5_file['variants']
        variant_data = {}
        
        for dataset_name in variants_group.keys():
            dataset = variants_group[dataset_name]
            if dataset.dtype.kind == 'S':  # String data
                variant_data[dataset_name] = [s.decode('utf-8') for s in dataset[:]]
            else:
                variant_data[dataset_name] = dataset[:]
        
        return pd.DataFrame(variant_data)
    
    def get_samples_df(self) -> pd.DataFrame:
        """
        Get sample information as a pandas DataFrame.
        
        Returns:
            DataFrame with sample information (columns from .fam file)
        """
        if self._h5_file is None:
            raise RuntimeError("HDF5 file is not open")
        
        if 'samples' not in self._h5_file:
            raise ValueError("No sample information found in HDF5 file")
        
        samples_group = self._h5_file['samples']
        sample_data = {}
        
        for dataset_name in samples_group.keys():
            dataset = samples_group[dataset_name]
            if dataset.dtype.kind == 'S':  # String data
                sample_data[dataset_name] = [s.decode('utf-8') for s in dataset[:]]
            else:
                sample_data[dataset_name] = dataset[:]
        
        return pd.DataFrame(sample_data)
    
    def get_variant_by_id(self, variant_id: str) -> Optional[int]:
        """
        Get the index of a variant by its ID.
        
        Args:
            variant_id: Variant ID to search for
            
        Returns:
            Index of the variant or None if not found
        """
        variants_df = self.get_variants_df()
        
        # Try to find in the 'snp' column (typical PLINK variant ID column)
        if 'snp' in variants_df.columns:
            matches = variants_df[variants_df['snp'] == variant_id].index
            if len(matches) > 0:
                return matches[0]
        
        # Try other possible ID columns
        for col in ['id', 'variant_id', 'rsid']:
            if col in variants_df.columns:
                matches = variants_df[variants_df[col] == variant_id].index
                if len(matches) > 0:
                    return matches[0]
        
        return None
    
    def get_sample_by_id(self, sample_id: str, id_column: str = 'iid') -> Optional[int]:
        """
        Get the index of a sample by its ID.
        
        Args:
            sample_id: Sample ID to search for
            id_column: Column name to search in (default: 'iid')
            
        Returns:
            Index of the sample or None if not found
        """
        samples_df = self.get_samples_df()
        
        if id_column not in samples_df.columns:
            raise ValueError(f"Column '{id_column}' not found in sample data")
        
        matches = samples_df[samples_df[id_column] == sample_id].index
        if len(matches) > 0:
            return matches[0]
        
        return None
    
    def subset_by_variant_ids(self, variant_ids: list) -> np.ndarray:
        """
        Get genotypes for specific variants by ID.
        
        Args:
            variant_ids: List of variant IDs
            
        Returns:
            Genotype matrix subset (n_samples, n_selected_variants)
        """
        indices = []
        for variant_id in variant_ids:
            idx = self.get_variant_by_id(variant_id)
            if idx is not None:
                indices.append(idx)
        
        if not indices:
            raise ValueError("No matching variants found")
        
        return self.get_genotypes(variant_indices=np.array(indices))
    
    def subset_by_sample_ids(self, sample_ids: list, id_column: str = 'iid') -> np.ndarray:
        """
        Get genotypes for specific samples by ID.
        
        Args:
            sample_ids: List of sample IDs
            id_column: Column name to search in
            
        Returns:
            Genotype matrix subset (n_selected_samples, n_variants)
        """
        indices = []
        for sample_id in sample_ids:
            idx = self.get_sample_by_id(sample_id, id_column)
            if idx is not None:
                indices.append(idx)
        
        if not indices:
            raise ValueError("No matching samples found")
        
        return self.get_genotypes(sample_indices=np.array(indices))
    
    def compute_allele_frequencies(self) -> np.ndarray:
        """
        Compute allele frequencies for all variants.
        
        Returns:
            Array of allele frequencies (n_variants,)
        """
        genotypes = self.get_genotypes()
        
        # Count non-missing genotypes (0 indicates missing in PLINK encoding)
        non_missing_mask = genotypes != 0
        
        # Calculate allele frequencies (genotype values: 1=het, 2=hom_alt)
        # Frequency = (number of alt alleles) / (2 * number of non-missing genotypes)
        alt_allele_counts = np.sum(genotypes - 1, axis=0)  # Convert 1,2 to 0,1 then sum
        total_alleles = 2 * np.sum(non_missing_mask, axis=0)
        
        # Avoid division by zero
        frequencies = np.divide(alt_allele_counts, total_alleles, 
                              out=np.zeros_like(alt_allele_counts, dtype=float), 
                              where=total_alleles != 0)
        
        return frequencies
    
    def info(self) -> str:
        """
        Get a summary of the HDF5 file contents.
        
        Returns:
            String with file information
        """
        metadata = self.metadata
        
        info_str = f"HDF5 PLINK File: {self.h5_file_path.name}\n"
        info_str += f"{'='*50}\n"
        info_str += f"Samples: {self.n_samples:,}\n"
        info_str += f"Variants: {self.n_variants:,}\n"
        
        if 'source_bed_file' in metadata:
            info_str += f"Source BED file: {metadata['source_bed_file']}\n"
        
        if 'conversion_tool' in metadata:
            info_str += f"Conversion tool: {metadata['conversion_tool']}\n"
        
        # Add encoding information
        if 'encoding' in metadata:
            info_str += f"Encoding: {metadata['encoding']}\n"
        
        return info_str


def read_h5_plink(h5_file_path: str) -> H5PLINKReader:
    """
    Convenience function to create an H5PLINKReader.
    
    Args:
        h5_file_path: Path to the HDF5 file
        
    Returns:
        H5PLINKReader instance
    """
    return H5PLINKReader(h5_file_path)


def list_h5_files(directory: str) -> list:
    """
    List all HDF5 files in a directory.
    
    Args:
        directory: Directory to search
        
    Returns:
        List of HDF5 file paths
    """
    directory_path = Path(directory)
    if not directory_path.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    return list(directory_path.rglob("*.h5"))
