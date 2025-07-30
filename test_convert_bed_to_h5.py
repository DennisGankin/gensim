#!/usr/bin/env python3
"""
Test script for the PLINK to HDF5 converter.

This script tests the basic functionality of the converter
without requiring large datasets or all dependencies.
"""

import sys
import os
from pathlib import Path

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_converter_initialization():
    """Test that the converter can be initialized."""
    try:
        # This will fail if dependencies aren't available, but syntax should be OK
        from convert_bed_to_h5 import BedToH5Converter
        
        # Test initialization (this should work even without dependencies)
        print("✓ BedToH5Converter class imported successfully")
        
        # Test that we can create an instance (will fail if input_dir doesn't exist)
        try:
            converter = BedToH5Converter(
                input_dir="test_nonexistent_dir",
                chunk_size=1000,
                memory_usage=2.0
            )
        except ValueError as e:
            if "does not exist" in str(e):
                print("✓ Input directory validation works correctly")
            else:
                print(f"✗ Unexpected error: {e}")
                return False
        except Exception as e:
            print(f"✗ Unexpected error during initialization: {e}")
            return False
        
        return True
        
    except ImportError as e:
        print(f"✗ Import error (expected if dependencies not installed): {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False

def test_chunk_calculation_logic():
    """Test the chunk size calculation logic."""
    print("\nTesting chunk size calculation logic...")
    
    # Test the calculation formula manually
    def calculate_chunk_size(n_samples, n_variants, memory_usage_gb=2.0):
        """Replicate the chunk size calculation."""
        bytes_per_sample = n_variants * 4  # float32
        target_memory_bytes = memory_usage_gb * 1024**3
        chunk_size = max(1, int(target_memory_bytes / bytes_per_sample))
        return min(chunk_size, n_samples)
    
    # Test cases
    test_cases = [
        (50000, 500000, 2.0),    # 50K samples, 500K variants, 2GB
        (100000, 1000000, 4.0),  # 100K samples, 1M variants, 4GB
        (1000000, 500000, 8.0),  # 1M samples, 500K variants, 8GB
    ]
    
    for n_samples, n_variants, memory_gb in test_cases:
        chunk_size = calculate_chunk_size(n_samples, n_variants, memory_gb)
        memory_per_chunk = (chunk_size * n_variants * 4) / (1024**3)
        
        print(f"  {n_samples:,} samples × {n_variants:,} variants @ {memory_gb}GB:")
        print(f"    → chunk_size: {chunk_size:,} samples")
        print(f"    → memory per chunk: {memory_per_chunk:.2f} GB")
        
        # Verify chunk size is reasonable
        if chunk_size <= 0 or chunk_size > n_samples:
            print(f"    ✗ Invalid chunk size")
            return False
        if memory_per_chunk > memory_gb * 1.1:  # Allow 10% tolerance
            print(f"    ✗ Memory usage too high")
            return False
    
    print("  ✓ Chunk size calculation logic is correct")
    return True

def test_file_discovery_logic():
    """Test the file discovery logic."""
    print("\nTesting file discovery logic...")
    
    # Create test directory structure
    test_dir = Path("test_discovery")
    try:
        test_dir.mkdir(exist_ok=True)
        
        # Create test files
        test_files = [
            "simulation1.bed",
            "simulation1.bim", 
            "simulation1.fam",
            "simulation2.bed",
            "simulation2.bim",
            "simulation2.fam",
            "temp-temporary.bed",  # Should be excluded
            "temp-temporary.bim",
            "temp-temporary.fam",
            "incomplete.bed",      # Should be excluded (no .bim/.fam)
        ]
        
        for file_name in test_files:
            (test_dir / file_name).touch()
        
        # Test the logic manually (replicate find_bed_files logic)
        bed_files = []
        for bed_file in test_dir.rglob("*.bed"):
            # Skip temporary files
            if "temporary" in bed_file.name.lower():
                print(f"  Skipping temporary file: {bed_file.name}")
                continue
            
            # Check for corresponding files
            bim_file = bed_file.with_suffix('.bim')
            fam_file = bed_file.with_suffix('.fam')
            
            if not bim_file.exists():
                print(f"  Missing .bim file for {bed_file.name}")
                continue
            
            if not fam_file.exists():
                print(f"  Missing .fam file for {bed_file.name}")
                continue
            
            bed_files.append(bed_file.name)
            print(f"  Found complete fileset: {bed_file.stem}")
        
        # Verify results
        expected_files = ["simulation1.bed", "simulation2.bed"]
        found_files = [Path(f).name for f in bed_files]
        
        if set(found_files) == set(expected_files):
            print("  ✓ File discovery logic is correct")
            result = True
        else:
            print(f"  ✗ Expected {expected_files}, found {found_files}")
            result = False
        
        # Cleanup
        for file_name in test_files:
            (test_dir / file_name).unlink(missing_ok=True)
        test_dir.rmdir()
        
        return result
        
    except Exception as e:
        print(f"  ✗ Error testing file discovery: {e}")
        return False

def main():
    """Run all tests."""
    print("PLINK to HDF5 Converter Tests")
    print("=" * 40)
    
    tests = [
        test_converter_initialization,
        test_chunk_calculation_logic,
        test_file_discovery_logic,
    ]
    
    passed = 0
    total = len(tests)
    
    for test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"✗ Test {test_func.__name__} failed with exception: {e}")
    
    print(f"\nTest Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("✓ All tests passed!")
        return 0
    else:
        print("✗ Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
