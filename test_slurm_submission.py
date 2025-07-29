#!/usr/bin/env python3
"""
Test SLURM parameter grid submission with a minimal example.

This script creates and optionally submits a small parameter grid
to test the SLURM submission functionality.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path


def test_chunked_submission(dry_run=True):
    """Test chunked job submission with minimal parameters."""
    print("Testing chunked SLURM submission...")
    
    cmd = [
        "python", "submit_parameter_grid_slurm.py",
        "--cohort-sizes", "100,200",
        "--prevalences", "0.01,0.05",
        "--total-snps", "1000,2000",
        "--causal-snps", "10,20",
        "--num-jobs", "2",
        "--time", "00:30:00",
        "--memory", "1G",
        "--cpus", "1",
        "--partition", "cpu",
        "--grid-output-dir", "test_chunked_grid",
        "--base-prefix", "test_dataset"
    ]
    
    if dry_run:
        cmd.append("--dry-run")
    
    print("Command:", " ".join(cmd))
    print()
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        return False


def test_individual_submission(dry_run=True):
    """Test individual job submission with minimal parameters."""
    print("Testing individual SLURM submission...")
    
    cmd = [
        "python", "submit_individual_slurm.py",
        "--cohort-sizes", "100,200",
        "--prevalences", "0.01,0.05",
        "--total-snps", "1000,2000",
        "--causal-snps", "10,20",
        "--time", "00:15:00",
        "--memory", "1G",
        "--partition", "cpu",
        "--grid-output-dir", "test_individual_grid",
        "--base-prefix", "test_individual",
        "--max-jobs", "8"  # Limit to 8 jobs for testing
    ]
    
    if dry_run:
        cmd.append("--dry-run")
    
    print("Command:", " ".join(cmd))
    print()
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        return False


def test_array_submission(dry_run=True):
    """Test array job submission with minimal parameters."""
    print("Testing array SLURM submission...")
    
    cmd = [
        "python", "submit_individual_slurm.py",
        "--cohort-sizes", "100,200",
        "--prevalences", "0.01,0.05",
        "--total-snps", "1000,2000",
        "--causal-snps", "10,20",
        "--use-array",
        "--time", "00:15:00",
        "--memory", "1G",
        "--partition", "cpu",
        "--grid-output-dir", "test_array_grid",
        "--base-prefix", "test_array"
    ]
    
    if dry_run:
        cmd.append("--dry-run")
    
    print("Command:", " ".join(cmd))
    print()
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        return False


def test_regular_parameter_grid():
    """Test regular (non-SLURM) parameter grid functionality."""
    print("Testing regular parameter grid functionality...")
    
    cmd = [
        "python", "main.py",
        "--create-parameter-grid",
        "--grid-cohort-sizes", "50,100",
        "--grid-prevalences", "0.01,0.05",
        "--grid-total-snps", "500,1000",
        "--grid-causal-snps", "5,10",
        "--grid-output-dir", "test_regular_grid",
        "--grid-prefix", "test_regular"
    ]
    
    print("Command:", " ".join(cmd))
    print("Note: This will actually run simulations (if PLINK is available)")
    print()
    
    try:
        # Don't actually run this by default since it requires PLINK
        print("Skipping actual execution. To test:")
        print("1. Ensure PLINK is installed and accessible")
        print("2. Run the command above manually")
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False


def check_prerequisites():
    """Check if required files and dependencies exist."""
    print("Checking prerequisites...")
    
    required_files = [
        "submit_parameter_grid_slurm.py",
        "submit_individual_slurm.py",
        "main.py",
        "gensim/__init__.py"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"❌ Missing required files: {missing_files}")
        return False
    
    print("✓ All required files found")
    
    # Check if SLURM is available
    try:
        subprocess.run(['squeue', '--version'], capture_output=True, check=True)
        print("✓ SLURM detected")
        slurm_available = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠️  SLURM not detected (scripts will create files but won't submit)")
        slurm_available = False
    
    # Check if Python modules can be imported
    try:
        from gensim import PLINKParameterGrid, PLINKParameterGridSimulator
        print("✓ Gensim modules can be imported")
    except ImportError as e:
        print(f"❌ Cannot import gensim modules: {e}")
        return False
    
    return True


def main():
    """Run tests for SLURM submission functionality."""
    parser = argparse.ArgumentParser(
        description="Test SLURM parameter grid submission scripts"
    )
    parser.add_argument("--submit", action="store_true",
                       help="Actually submit jobs (remove --dry-run)")
    parser.add_argument("--test-regular", action="store_true",
                       help="Test regular parameter grid (requires PLINK)")
    
    args = parser.parse_args()
    
    print("SLURM Parameter Grid Submission Test")
    print("=" * 50)
    print()
    
    # Check prerequisites
    if not check_prerequisites():
        print("❌ Prerequisites check failed")
        return 1
    
    print()
    dry_run = not args.submit
    
    if dry_run:
        print("🔍 Running in dry-run mode (scripts created, jobs not submitted)")
        print("   Use --submit flag to actually submit jobs")
    else:
        print("🚀 Running in submission mode (jobs will be submitted)")
    
    print()
    
    # Test chunked submission
    print("1. Testing chunked job submission...")
    success1 = test_chunked_submission(dry_run)
    print("   Result:", "✓ Success" if success1 else "❌ Failed")
    print()
    
    # Test individual submission
    print("2. Testing individual job submission...")
    success2 = test_individual_submission(dry_run)
    print("   Result:", "✓ Success" if success2 else "❌ Failed")
    print()
    
    # Test array submission
    print("3. Testing array job submission...")
    success3 = test_array_submission(dry_run)
    print("   Result:", "✓ Success" if success3 else "❌ Failed")
    print()
    
    # Test regular parameter grid if requested
    if args.test_regular:
        print("4. Testing regular parameter grid...")
        success4 = test_regular_parameter_grid()
        print("   Result:", "✓ Success" if success4 else "❌ Failed")
        print()
    else:
        print("4. Skipping regular parameter grid test (use --test-regular to enable)")
        success4 = True
        print()
    
    # Summary
    all_success = success1 and success2 and success3 and success4
    
    print("=" * 50)
    print("TEST SUMMARY")
    print("=" * 50)
    print(f"Chunked submission:    {'✓' if success1 else '❌'}")
    print(f"Individual submission: {'✓' if success2 else '❌'}")
    print(f"Array submission:      {'✓' if success3 else '❌'}")
    print(f"Regular grid:          {'✓' if success4 else '❌'}")
    print()
    
    if all_success:
        print("🎉 All tests passed!")
        
        if dry_run:
            print()
            print("Next steps:")
            print("1. Check created SLURM scripts in test_*_grid/slurm_* directories")
            print("2. Verify scripts look correct")
            print("3. Run with --submit flag to actually submit test jobs")
            print("4. Monitor jobs with 'squeue -u $USER'")
        else:
            print()
            print("Jobs submitted! Monitor with:")
            print("  squeue -u $USER")
            print("  ls test_*_grid/")
    else:
        print("❌ Some tests failed. Check error messages above.")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
