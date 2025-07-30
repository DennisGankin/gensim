#!/usr/bin/env python3
"""
Debug script to check memory allocation and dataset sizes.
"""

import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))

from convert_bed_to_h5 import BedToH5Converter
from pysnptools.snpreader import Bed

def main():
    print("=== SLURM Memory Debug ===")
    
    # Check SLURM environment
    slurm_vars = ['SLURM_JOB_ID', 'SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU', 'SLURM_CPUS_PER_TASK']
    for var in slurm_vars:
        value = os.environ.get(var, 'NOT SET')
        print(f"{var}: {value}")
    
    print("\n=== Memory Detection Test ===")
    converter = BedToH5Converter(input_dir="data/ukb_sim", output_dir="test_output")
    
    # Test memory detection
    slurm_memory = converter._get_slurm_memory_limit()
    print(f"Detected SLURM memory: {slurm_memory} GB")
    
    print("\n=== Dataset Analysis ===")
    bed_files = converter.find_bed_files()
    
    for i, (bed_file_path, prefix) in enumerate(bed_files[:3]):  # Check first 3 files
        print(f"\nFile {i+1}: {bed_file_path}")
        
        try:
            plink_prefix = str(bed_file_path.with_suffix(''))
            bed_reader = Bed(plink_prefix, count_A1=False)
            
            n_samples = bed_reader.iid_count
            n_variants = bed_reader.sid_count
            
            # Calculate memory requirement
            memory_gb = (n_samples * n_variants * 4) / (1024**3)
            
            print(f"  Samples: {n_samples:,}")
            print(f"  Variants: {n_variants:,}")
            print(f"  Memory needed: {memory_gb:.1f} GB")
            
            # Test chunk size calculation
            chunk_size = converter.calculate_chunk_size(n_samples, n_variants)
            print(f"  Calculated chunk size: {chunk_size:,}")
            
            # Test memory check
            needs_chunking, required_gb, available_gb = converter.check_memory_requirements(n_samples, n_variants)
            print(f"  Needs chunking: {needs_chunking}")
            print(f"  Available memory: {available_gb:.1f} GB")
            
        except Exception as e:
            print(f"  ERROR: {e}")

if __name__ == "__main__":
    main()
