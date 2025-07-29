#!/usr/bin/env python3
"""
SLURM job submission script for parallel parameter grid simulation.

This script splits a parameter grid into chunks and submits SLURM jobs
to run each chunk in parallel across the cluster.
"""

import os
import sys
import argparse
import json
import subprocess
import math
from pathlib import Path
from typing import List, Dict, Any, Tuple

from gensim import PLINKParameterGrid


def create_slurm_script(
    job_name: str,
    output_file: str,
    error_file: str,
    time_limit: str = "02:00:00",
    memory: str = "4G",
    cpus: int = 1,
    partition: str = "cpu",
    account: str = None,
    commands: List[str] = None
) -> str:
    """
    Create a SLURM job script.
    
    Args:
        job_name: Name of the SLURM job
        output_file: Path for stdout output
        error_file: Path for stderr output
        time_limit: Job time limit (HH:MM:SS format)
        memory: Memory requirement
        cpus: Number of CPUs
        partition: SLURM partition
        account: SLURM account (optional)
        commands: List of commands to execute
        
    Returns:
        SLURM script content as string
    """
    script_lines = [
        "#!/bin/bash",
        "",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={output_file}",
        f"#SBATCH --error={error_file}",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --mem-per-cpu={memory}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --partition={partition}",
    ]
    
    if account:
        script_lines.append(f"#SBATCH --account={account}")
    
    script_lines.extend([
        "",
        "# Load modules and activate mamba environment",
        "# module load mamba  # Uncomment if needed on your cluster",
        "source ~/.bashrc",
        "mamba activate gensim",
        "",
        "# Set working directory",
        f"cd {os.getcwd()}",
        "",
        "# Print job information",
        "echo \"Job started at: $(date)\"",
        "echo \"Job ID: $SLURM_JOB_ID\"",
        "echo \"Node: $SLURM_NODELIST\"",
        "echo \"Working directory: $(pwd)\"",
        "",
        "# Run commands",
    ])
    
    if commands:
        script_lines.extend(commands)
    
    script_lines.extend([
        "",
        "echo \"Job finished at: $(date)\"",
    ])
    
    return "\n".join(script_lines)


def split_parameter_combinations(
    grid_config: PLINKParameterGrid,
    num_chunks: int
) -> List[Tuple[int, int, List[Dict[str, Any]]]]:
    """
    Split parameter combinations into chunks for parallel processing.
    
    Args:
        grid_config: Parameter grid configuration
        num_chunks: Number of chunks to create
        
    Returns:
        List of tuples: (start_idx, end_idx, combinations_chunk)
    """
    combinations = grid_config.get_parameter_combinations()
    total_combinations = len(combinations)
    
    if num_chunks > total_combinations:
        num_chunks = total_combinations
        print(f"Warning: More chunks ({num_chunks}) than combinations ({total_combinations}). "
              f"Reducing to {total_combinations} chunks.")
    
    chunk_size = math.ceil(total_combinations / num_chunks)
    chunks = []
    
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, total_combinations)
        
        if start_idx >= total_combinations:
            break
            
        chunk_combinations = combinations[start_idx:end_idx]
        chunks.append((start_idx, end_idx - 1, chunk_combinations))
    
    return chunks


def create_chunk_parameter_file(
    chunk_combinations: List[Dict[str, Any]],
    chunk_id: int,
    base_grid_config: PLINKParameterGrid,
    output_dir: str
) -> str:
    """
    Create a parameter file for a specific chunk.
    
    Args:
        chunk_combinations: Parameter combinations for this chunk
        chunk_id: Chunk identifier
        base_grid_config: Base grid configuration
        output_dir: Output directory for chunk files
        
    Returns:
        Path to created parameter file
    """
    # Extract unique values for each parameter from the chunk
    cohort_sizes = sorted(list(set(combo['cohort_size'] for combo in chunk_combinations)))
    prevalences = sorted(list(set(combo['prevalence'] for combo in chunk_combinations)))
    total_snps = sorted(list(set(combo['total_snps'] for combo in chunk_combinations)))
    causal_snps = sorted(list(set(combo['causal_snps'] for combo in chunk_combinations)))
    
    chunk_config = {
        'cohort_sizes': cohort_sizes,
        'prevalences': prevalences,
        'total_snps': total_snps,
        'causal_snps': causal_snps,
        'grid_output_dir': base_grid_config.grid_output_dir,  # No chunk subdirectory
        'base_prefix': base_grid_config.base_prefix,  # Keep original prefix
        'case_control_ratio': base_grid_config.case_control_ratio,
        'min_freq': base_grid_config.min_freq,
        'max_freq': base_grid_config.max_freq,
        'het_odds_ratio': base_grid_config.het_odds_ratio,
        'hom_odds_ratio': base_grid_config.hom_odds_ratio,
        'plink_executable': base_grid_config.plink_executable,
        'random_seed': base_grid_config.random_seed
    }
    
    chunk_param_file = os.path.join(output_dir, f"chunk_{chunk_id}_params.json")
    
    with open(chunk_param_file, 'w') as f:
        json.dump(chunk_config, f, indent=2)
    
    return chunk_param_file


def submit_slurm_jobs(
    grid_config: PLINKParameterGrid,
    num_jobs: int,
    slurm_args: Dict[str, Any],
    dry_run: bool = False
) -> List[str]:
    """
    Submit SLURM jobs for parameter grid simulation.
    
    Args:
        grid_config: Parameter grid configuration
        num_jobs: Number of parallel jobs to submit
        slurm_args: SLURM job arguments
        dry_run: If True, create scripts but don't submit jobs
        
    Returns:
        List of job IDs (empty if dry_run=True)
    """
    # Create output directories
    slurm_output_dir = os.path.join(grid_config.grid_output_dir, "slurm_jobs")
    os.makedirs(slurm_output_dir, exist_ok=True)
    
    # Split combinations into chunks
    chunks = split_parameter_combinations(grid_config, num_jobs)
    
    print(f"Splitting {len(grid_config.get_parameter_combinations())} combinations into {len(chunks)} chunks:")
    for i, (start_idx, end_idx, chunk_combos) in enumerate(chunks):
        print(f"  Chunk {i+1}: combinations {start_idx+1}-{end_idx+1} ({len(chunk_combos)} combinations)")
    
    job_ids = []
    
    for chunk_id, (start_idx, end_idx, chunk_combinations) in enumerate(chunks):
        # Create parameter file for this chunk
        chunk_param_file = create_chunk_parameter_file(
            chunk_combinations, chunk_id + 1, grid_config, slurm_output_dir
        )
        
        # Create SLURM script
        job_name = f"plink_grid_chunk_{chunk_id + 1}"
        output_file = os.path.join(slurm_output_dir, f"{job_name}.out")
        error_file = os.path.join(slurm_output_dir, f"{job_name}.err")
        
        # Build command to run the chunk
        chunk_commands = [
            f"echo \"Processing chunk {chunk_id + 1} with {len(chunk_combinations)} combinations\"",
            "",
            "# Run parameter grid for this chunk",
            f"python -c \"",
            "import json",
            "from gensim import PLINKParameterGrid, PLINKParameterGridSimulator",
            f"with open('{chunk_param_file}', 'r') as f:",
            "    config_dict = json.load(f)",
            "config_dict['hom_odds_ratio'] = config_dict['hom_odds_ratio'] if config_dict['hom_odds_ratio'] == 'mult' else float(config_dict['hom_odds_ratio'])",
            "grid_config = PLINKParameterGrid(**config_dict)",
            "simulator = PLINKParameterGridSimulator(grid_config)",
            "results = simulator.run_parameter_grid()",
            "print(f'Chunk results: {results[\\\"successful\\\"]}/{results[\\\"total_combinations\\\"]} successful')",
            "\""
        ]
        
        slurm_script_content = create_slurm_script(
            job_name=job_name,
            output_file=output_file,
            error_file=error_file,
            time_limit=slurm_args.get('time', '02:00:00'),
            memory=slurm_args.get('memory', '4G'),
            cpus=slurm_args.get('cpus', 1),
            partition=slurm_args.get('partition', 'cpu'),
            account=slurm_args.get('account'),
            commands=chunk_commands
        )
        
        # Write SLURM script
        script_file = os.path.join(slurm_output_dir, f"{job_name}.sh")
        with open(script_file, 'w') as f:
            f.write(slurm_script_content)
        
        # Make script executable
        os.chmod(script_file, 0o755)
        
        if dry_run:
            print(f"Created SLURM script: {script_file}")
        else:
            # Submit job
            try:
                result = subprocess.run(
                    ['sbatch', script_file],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Extract job ID from sbatch output
                job_id = result.stdout.strip().split()[-1]
                job_ids.append(job_id)
                
                print(f"Submitted job {job_id}: {job_name}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error submitting job {job_name}: {e}")
                print(f"STDERR: {e.stderr}")
    
    return job_ids


def create_monitor_script(grid_output_dir: str, job_ids: List[str]) -> str:
    """
    Create a monitoring script to check job status and collect results.
    
    Args:
        grid_output_dir: Grid output directory
        job_ids: List of submitted job IDs
        
    Returns:
        Path to monitoring script
    """
    monitor_script = os.path.join(grid_output_dir, "monitor_jobs.sh")
    
    script_content = f"""#!/bin/bash

# Monitor SLURM jobs for parameter grid simulation
JOB_IDS=({' '.join(job_ids)})
GRID_DIR="{grid_output_dir}"

echo "Monitoring {len(job_ids)} SLURM jobs..."
echo "Job IDs: {' '.join(job_ids)}"
echo ""

# Function to check job status
check_jobs() {{
    echo "Checking job status at $(date):"
    for job_id in "${{JOB_IDS[@]}}"; do
        status=$(squeue -j $job_id -h -o "%T" 2>/dev/null || echo "COMPLETED/FAILED")
        echo "  Job $job_id: $status"
    done
    echo ""
}}

# Function to collect results
collect_results() {{
    echo "Collecting results from SLURM job outputs..."
    
    total_combinations=0
    total_successful=0
    total_failed=0
    
    # Check each SLURM job output file for results
    for output_file in "$GRID_DIR"/slurm_jobs/plink_grid_chunk_*.out; do
        if [ -f "$output_file" ]; then
            echo "Checking results in: $(basename $output_file)"
            
            # Count successful and failed combinations from output
            successful_count=$(grep -c "✓ Simulation completed successfully!" "$output_file" 2>/dev/null || echo "0")
            failed_count=$(grep -c "✗ Simulation failed!" "$output_file" 2>/dev/null || echo "0")
            chunk_total=$((successful_count + failed_count))
            
            if [ $chunk_total -gt 0 ]; then
                echo "  Combinations: $chunk_total, Successful: $successful_count, Failed: $failed_count"
                total_combinations=$((total_combinations + chunk_total))
                total_successful=$((total_successful + successful_count))
                total_failed=$((total_failed + failed_count))
            fi
        fi
    done
    
    # Also count actual dataset directories created
    dataset_count=$(find "$GRID_DIR" -maxdepth 1 -name "dataset_*" -type d | wc -l)
    
    echo ""
    echo "OVERALL RESULTS:"
    echo "Total combinations processed: $total_combinations"
    echo "Successful simulations: $total_successful"
    echo "Failed simulations: $total_failed"
    echo "Dataset directories created: $dataset_count"
    if [ $total_combinations -gt 0 ]; then
        echo "Success rate: $(echo "scale=1; $total_successful * 100 / $total_combinations" | bc -l 2>/dev/null || echo "N/A")%"
    fi
    echo ""
    echo "Generated datasets are in: $GRID_DIR/dataset_*/"
    echo "SLURM job outputs in: $GRID_DIR/slurm_jobs/"
}}

# Check initial status
check_jobs

# Wait for jobs to complete
echo "Waiting for jobs to complete..."
echo "You can also check status manually with: squeue -u \\$USER"
echo ""

while true; do
    running_jobs=0
    for job_id in "${{JOB_IDS[@]}}"; do
        if squeue -j $job_id &>/dev/null; then
            running_jobs=$((running_jobs + 1))
        fi
    done
    
    if [ $running_jobs -eq 0 ]; then
        echo "All jobs completed!"
        check_jobs
        collect_results
        break
    else
        echo "Still running: $running_jobs jobs"
        sleep 60  # Check every minute
    fi
done
"""
    
    with open(monitor_script, 'w') as f:
        f.write(script_content)
    
    os.chmod(monitor_script, 0o755)
    return monitor_script


def main():
    """Main function for SLURM job submission."""
    parser = argparse.ArgumentParser(
        description="Submit SLURM jobs for parallel parameter grid simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit 8 parallel jobs with basic parameters
  python submit_parameter_grid_slurm.py \\
    --cohort-sizes "500,1000,2000" \\
    --prevalences "0.01,0.05,0.1" \\
    --total-snps "5000,10000" \\
    --causal-snps "50,100" \\
    --num-jobs 8
  
  # Submit with custom SLURM parameters
  python submit_parameter_grid_slurm.py \\
    --cohort-sizes "1000,2000" \\
    --prevalences "0.01,0.1" \\
    --total-snps "10000" \\
    --causal-snps "100,200" \\
    --num-jobs 4 \\
    --time "04:00:00" \\
    --memory "8G" \\
    --partition "gpu" \\
    --account "myaccount"
  
  # Dry run to create scripts without submitting
  python submit_parameter_grid_slurm.py \\
    --cohort-sizes "500,1000" \\
    --prevalences "0.01" \\
    --total-snps "5000" \\
    --causal-snps "50" \\
    --num-jobs 2 \\
    --dry-run
        """
    )
    
    # Grid parameters
    parser.add_argument("--cohort-sizes", required=True,
                       help="Comma-separated list of cohort sizes")
    parser.add_argument("--prevalences", required=True,
                       help="Comma-separated list of prevalences")
    parser.add_argument("--total-snps", required=True,
                       help="Comma-separated list of total SNP counts")
    parser.add_argument("--causal-snps", required=True,
                       help="Comma-separated list of causal SNP counts")
    parser.add_argument("--grid-output-dir", default="parameter_grid_slurm",
                       help="Output directory for parameter grid")
    parser.add_argument("--base-prefix", default="dataset",
                       help="Base prefix for dataset names")
    
    # Grid configuration
    parser.add_argument("--case-control-ratio", type=float, default=1.0,
                       help="Case-to-control ratio")
    parser.add_argument("--min-freq", type=float, default=0.01,
                       help="Minimum allele frequency")
    parser.add_argument("--max-freq", type=float, default=0.5,
                       help="Maximum allele frequency")
    parser.add_argument("--het-or", type=float, default=1.5,
                       help="Heterozygote odds ratio for causal SNPs")
    parser.add_argument("--hom-or", default="mult",
                       help="Homozygote odds ratio for causal SNPs")
    parser.add_argument("--plink-executable", default="plink",
                       help="PLINK executable name or path")
    parser.add_argument("--random-seed", type=int,
                       help="Random seed for reproducibility")
    
    # SLURM parameters
    parser.add_argument("--num-jobs", type=int, required=True,
                       help="Number of parallel SLURM jobs to submit")
    parser.add_argument("--time", default="02:00:00",
                       help="SLURM time limit (HH:MM:SS)")
    parser.add_argument("--memory", default="4G",
                       help="SLURM memory requirement")
    parser.add_argument("--cpus", type=int, default=1,
                       help="Number of CPUs per job")
    parser.add_argument("--partition", default="cpu",
                       help="SLURM partition")
    parser.add_argument("--account",
                       help="SLURM account")
    
    # Options
    parser.add_argument("--dry-run", action="store_true",
                       help="Create scripts but don't submit jobs")
    
    args = parser.parse_args()
    
    # Parse parameter lists
    def parse_int_list(value_str):
        return [int(x.strip()) for x in value_str.split(',')]
    
    def parse_float_list(value_str):
        return [float(x.strip()) for x in value_str.split(',')]
    
    # Parse homozygote odds ratio
    hom_or = args.hom_or
    if hom_or != "mult":
        try:
            hom_or = float(hom_or)
        except ValueError:
            print(f"Error: Invalid homozygote odds ratio: {hom_or}. Must be 'mult' or a number.")
            return 1
    
    # Create parameter grid configuration
    try:
        grid_config = PLINKParameterGrid(
            cohort_sizes=parse_int_list(args.cohort_sizes),
            prevalences=parse_float_list(args.prevalences),
            total_snps=parse_int_list(args.total_snps),
            causal_snps=parse_int_list(args.causal_snps),
            grid_output_dir=args.grid_output_dir,
            base_prefix=args.base_prefix,
            case_control_ratio=args.case_control_ratio,
            min_freq=args.min_freq,
            max_freq=args.max_freq,
            het_odds_ratio=args.het_or,
            hom_odds_ratio=hom_or,
            plink_executable=args.plink_executable,
            random_seed=args.random_seed
        )
    except Exception as e:
        print(f"Error creating parameter grid configuration: {e}")
        return 1
    
    # Display grid information
    total_combinations = len(grid_config.get_parameter_combinations())
    print(f"Parameter Grid Configuration:")
    print(f"  Cohort sizes: {grid_config.cohort_sizes}")
    print(f"  Prevalences: {grid_config.prevalences}")
    print(f"  Total SNPs: {grid_config.total_snps}")
    print(f"  Causal SNPs: {grid_config.causal_snps}")
    print(f"  Total combinations: {total_combinations}")
    print(f"  Output directory: {grid_config.grid_output_dir}")
    print()
    
    # SLURM configuration
    slurm_args = {
        'time': args.time,
        'memory': args.memory,
        'cpus': args.cpus,
        'partition': args.partition,
        'account': args.account
    }
    
    print(f"SLURM Configuration:")
    print(f"  Number of jobs: {args.num_jobs}")
    print(f"  Time limit: {args.time}")
    print(f"  Memory: {args.memory}")
    print(f"  CPUs per job: {args.cpus}")
    print(f"  Partition: {args.partition}")
    if args.account:
        print(f"  Account: {args.account}")
    print(f"  Dry run: {args.dry_run}")
    print()
    
    # Submit jobs
    try:
        job_ids = submit_slurm_jobs(grid_config, args.num_jobs, slurm_args, args.dry_run)
        
        if not args.dry_run and job_ids:
            print(f"\\nSuccessfully submitted {len(job_ids)} jobs!")
            print(f"Job IDs: {', '.join(job_ids)}")
            
            # Create monitoring script
            monitor_script = create_monitor_script(args.grid_output_dir, job_ids)
            print(f"\\nMonitoring script created: {monitor_script}")
            print(f"To monitor jobs, run: {monitor_script}")
            print(f"To check job status manually: squeue -u $USER")
            
        elif args.dry_run:
            print(f"\\nDry run completed. SLURM scripts created in:")
            print(f"  {os.path.join(args.grid_output_dir, 'slurm_jobs')}")
            print(f"To submit jobs, remove --dry-run flag and run again.")
            
        return 0
        
    except Exception as e:
        print(f"Error submitting SLURM jobs: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
