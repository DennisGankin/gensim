#!/usr/bin/env python3
"""
Examples of SLURM-based parameter grid simulation.

This script demonstrates different ways to submit parameter grid simulations
to SLURM clusters for parallel execution.
"""

import os
import subprocess
import sys
from pathlib import Path


def example_chunked_jobs():
    """Example of submitting chunked parameter grid jobs."""
    print("=" * 60)
    print("CHUNKED SLURM JOBS EXAMPLE")
    print("=" * 60)
    
    print("This approach splits the parameter grid into chunks and submits")
    print("a separate SLURM job for each chunk. Good for moderate-sized grids.")
    print()
    
    # Example command
    cmd = [
        "python", "submit_parameter_grid_slurm.py",
        "--cohort-sizes", "500,1000,2000",
        "--prevalences", "0.01,0.05,0.1",
        "--total-snps", "5000,10000,20000",
        "--causal-snps", "50,100,200",
        "--num-jobs", "8",
        "--time", "02:00:00",
        "--memory", "4G",
        "--partition", "cpu",
        "--grid-output-dir", "chunked_grid_example",
        "--dry-run"  # Remove this to actually submit
    ]
    
    print("Example command:")
    print(" ".join(cmd))
    print()
    
    print("This creates a 3×3×3×3 = 81 combination grid split into 8 chunks")
    print("Each chunk will process ~10-11 combinations")
    print()
    
    # Show what this would do
    grid_size = 3 * 3 * 3 * 3
    num_jobs = 8
    combinations_per_job = grid_size // num_jobs
    
    print(f"Grid breakdown:")
    print(f"  Total combinations: {grid_size}")
    print(f"  Number of SLURM jobs: {num_jobs}")
    print(f"  Combinations per job: ~{combinations_per_job}")
    print(f"  Expected runtime: ~2 hours per job")
    print()


def example_individual_jobs():
    """Example of submitting individual jobs for each combination."""
    print("=" * 60)
    print("INDIVIDUAL SLURM JOBS EXAMPLE")
    print("=" * 60)
    
    print("This approach submits one SLURM job per parameter combination.")
    print("Maximum parallelization but creates many jobs. Good for large grids.")
    print()
    
    # Example command for smaller grid
    cmd = [
        "python", "submit_individual_slurm.py",
        "--cohort-sizes", "1000,2000",
        "--prevalences", "0.01,0.1",
        "--total-snps", "10000,20000",
        "--causal-snps", "100,200",
        "--time", "30:00",
        "--memory", "2G",
        "--partition", "cpu",
        "--grid-output-dir", "individual_jobs_example",
        "--dry-run"
    ]
    
    print("Example command:")
    print(" ".join(cmd))
    print()
    
    grid_size = 2 * 2 * 2 * 2
    print(f"This creates {grid_size} individual SLURM jobs")
    print("Each job processes exactly 1 combination")
    print("Jobs can run simultaneously (limited by cluster capacity)")
    print()


def example_array_jobs():
    """Example of using SLURM array jobs."""
    print("=" * 60)
    print("SLURM ARRAY JOBS EXAMPLE")
    print("=" * 60)
    
    print("This approach uses SLURM array jobs for efficient job submission.")
    print("Single submission creates many parallel tasks. Most efficient for very large grids.")
    print()
    
    cmd = [
        "python", "submit_individual_slurm.py",
        "--cohort-sizes", "500,1000,2000,5000",
        "--prevalences", "0.001,0.01,0.05,0.1,0.2",
        "--total-snps", "5000,10000,20000,50000",
        "--causal-snps", "50,100,200,500",
        "--use-array",
        "--time", "01:00:00",
        "--memory", "3G",
        "--partition", "cpu",
        "--grid-output-dir", "array_jobs_example",
        "--dry-run"
    ]
    
    print("Example command:")
    print(" ".join(cmd))
    print()
    
    grid_size = 4 * 5 * 4 * 4
    print(f"This creates 1 array job with {grid_size} tasks")
    print("SLURM manages task scheduling automatically")
    print("More efficient than individual job submission")
    print()


def example_custom_slurm_parameters():
    """Example with custom SLURM parameters for specific cluster setup."""
    print("=" * 60)
    print("CUSTOM SLURM PARAMETERS EXAMPLE")
    print("=" * 60)
    
    print("This shows how to customize SLURM parameters for your specific cluster.")
    print()
    
    cmd = [
        "python", "submit_parameter_grid_slurm.py",
        "--cohort-sizes", "1000,5000,10000",
        "--prevalences", "0.01,0.1",
        "--total-snps", "50000,100000",
        "--causal-snps", "500,1000",
        "--num-jobs", "4",
        # Custom SLURM parameters
        "--time", "08:00:00",        # 8 hours
        "--memory", "16G",           # 16GB RAM
        "--cpus", "4",               # 4 CPUs per job
        "--partition", "highmem",    # High-memory partition
        "--account", "myproject",    # Billing account
        "--grid-output-dir", "large_grid_example",
        "--dry-run"
    ]
    
    print("Example command for large-scale simulation:")
    print(" ".join(cmd))
    print()
    
    print("Custom parameters explained:")
    print("  --time 08:00:00     : 8-hour time limit for large datasets")
    print("  --memory 16G        : 16GB RAM for memory-intensive simulations")
    print("  --cpus 4            : 4 CPUs per job (if PLINK supports parallelization)")
    print("  --partition highmem : Use high-memory partition")
    print("  --account myproject : Charge compute time to specific account")
    print()


def example_monitoring_jobs():
    """Example of monitoring submitted jobs."""
    print("=" * 60)
    print("JOB MONITORING EXAMPLE")
    print("=" * 60)
    
    print("After submitting jobs, you can monitor them in several ways:")
    print()
    
    print("1. Automatic monitoring script (created by chunked submission):")
    print("   ./parameter_grid_output/monitor_jobs.sh")
    print()
    
    print("2. Manual SLURM commands:")
    print("   squeue -u $USER                    # Show your jobs")
    print("   squeue -j 12345                   # Show specific job")
    print("   sacct -j 12345                    # Show job accounting info")
    print("   scancel 12345                     # Cancel job")
    print("   scancel -u $USER                  # Cancel all your jobs")
    print()
    
    print("3. Check output files:")
    print("   ls parameter_grid_output/slurm_jobs/*.out    # Stdout files")
    print("   ls parameter_grid_output/slurm_jobs/*.err    # Stderr files")
    print()
    
    print("4. Collect results:")
    print("   find parameter_grid_output -name '*.bed' | wc -l    # Count generated datasets")
    print("   find parameter_grid_output -name 'grid_summary.txt' # Find summary reports")
    print()


def example_failure_recovery():
    """Example of handling and recovering from job failures."""
    print("=" * 60)
    print("FAILURE RECOVERY EXAMPLE")
    print("=" * 60)
    
    print("Sometimes jobs fail due to various reasons. Here's how to handle failures:")
    print()
    
    print("1. Identify failed jobs:")
    print("   sacct -u $USER --state=FAILED     # Show failed jobs")
    print("   grep -r 'FAILED' parameter_grid_output/slurm_jobs/*.err")
    print()
    
    print("2. Check failure reasons:")
    print("   # Common issues:")
    print("   # - Time limit exceeded (increase --time)")
    print("   # - Memory limit exceeded (increase --memory)")
    print("   # - PLINK executable not found (check --plink-executable)")
    print("   # - Disk space issues")
    print("   # - Network/filesystem problems")
    print()
    
    print("3. Resubmit failed combinations:")
    print("   # Option A: Resubmit entire failed chunk")
    print("   # Option B: Create new grid with only failed parameter combinations")
    print("   # Option C: Submit individual jobs for failed combinations")
    print()
    
    print("4. Example resubmission:")
    resubmit_cmd = [
        "python", "submit_individual_slurm.py",
        "--cohort-sizes", "5000",      # Only the failed combination
        "--prevalences", "0.2",
        "--total-snps", "100000",
        "--causal-snps", "1000",
        "--time", "12:00:00",          # Increased time limit
        "--memory", "32G",             # Increased memory
        "--grid-output-dir", "recovery_jobs"
    ]
    print("   " + " ".join(resubmit_cmd))
    print()


def show_best_practices():
    """Show best practices for SLURM parameter grid submission."""
    print("=" * 60)
    print("BEST PRACTICES")
    print("=" * 60)
    
    print("1. Start small:")
    print("   - Test with a small grid (2×2×2×2 = 16 combinations)")
    print("   - Use --dry-run flag to check scripts before submission")
    print("   - Verify one combination works before scaling up")
    print()
    
    print("2. Choose appropriate job size:")
    print("   - Small grids (< 50 combinations): Individual jobs")
    print("   - Medium grids (50-500 combinations): Chunked jobs (5-20 chunks)")
    print("   - Large grids (> 500 combinations): Array jobs")
    print()
    
    print("3. Resource allocation:")
    print("   - Start with conservative estimates (1-2 hours, 2-4GB RAM)")
    print("   - Monitor actual usage and adjust")
    print("   - Consider SNP count and sample size when estimating resources")
    print()
    
    print("4. Cluster etiquette:")
    print("   - Don't submit thousands of jobs simultaneously")
    print("   - Use appropriate partitions (don't use GPU partition for CPU-only work)")
    print("   - Set reasonable time limits (don't ask for 24 hours if you need 1 hour)")
    print("   - Clean up failed/cancelled jobs")
    print()
    
    print("5. Data management:")
    print("   - Organize output in clear directory structure")
    print("   - Keep track of parameter combinations and their outputs")
    print("   - Compress or archive completed datasets to save space")
    print("   - Use meaningful naming conventions")
    print()


def main():
    """Run all examples."""
    print("GENSIM SLURM Parameter Grid Examples")
    print("=" * 80)
    print("This script shows different approaches for running parameter grid")
    print("simulations on SLURM clusters.")
    print()
    
    # Check if SLURM is available
    try:
        subprocess.run(['squeue', '--version'], capture_output=True, check=True)
        slurm_available = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        slurm_available = False
        print("⚠️  SLURM not detected. Examples are for reference only.")
        print()
    
    example_chunked_jobs()
    print("\n")
    
    example_individual_jobs()
    print("\n")
    
    example_array_jobs()
    print("\n")
    
    example_custom_slurm_parameters()
    print("\n")
    
    example_monitoring_jobs()
    print("\n")
    
    example_failure_recovery()
    print("\n")
    
    show_best_practices()
    
    print("\n" + "=" * 80)
    if slurm_available:
        print("✓ SLURM detected. You can run these examples on your cluster.")
        print("Remove --dry-run flags to actually submit jobs.")
    else:
        print("ℹ️  Install SLURM or run on a SLURM cluster to use these examples.")
    
    print("\nFor more information:")
    print("  python submit_parameter_grid_slurm.py --help")
    print("  python submit_individual_slurm.py --help")


if __name__ == "__main__":
    main()
