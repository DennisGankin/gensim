#!/usr/bin/env python3
"""
Submit individual parameter combinations as separate SLURM jobs.

This script submits each parameter combination as a separate SLURM job,
providing maximum parallelization for large parameter grids.
"""

import os
import sys
import argparse
import subprocess
import json
from typing import List, Dict, Any

from gensim import PLINKParameterGrid


def create_single_combination_script(
    combination: Dict[str, Any],
    combination_id: int,
    base_config: PLINKParameterGrid,
    slurm_args: Dict[str, Any]
) -> str:
    """
    Create a SLURM script for a single parameter combination.
    
    Args:
        combination: Parameter combination dictionary
        combination_id: Unique ID for this combination
        base_config: Base parameter grid configuration
        slurm_args: SLURM job arguments
        
    Returns:
        Path to created SLURM script
    """
    # Create combination-specific output directory
    combo_name = base_config.get_combination_name(combination)
    combo_dir = os.path.join(base_config.grid_output_dir, combo_name)
    
    # SLURM script content
    job_name = f"plink_{combo_name}"
    slurm_output_dir = os.path.join(base_config.grid_output_dir, "slurm_logs")
    os.makedirs(slurm_output_dir, exist_ok=True)
    
    output_file = os.path.join(slurm_output_dir, f"{job_name}.out")
    error_file = os.path.join(slurm_output_dir, f"{job_name}.err")
    
    # Prepare values to avoid f-string issues
    hom_or_value = "'mult'" if base_config.hom_odds_ratio == 'mult' else str(base_config.hom_odds_ratio)
    random_seed_value = str(base_config.random_seed) if base_config.random_seed else 'None'
    
    script_content = f"""#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --output={output_file}
#SBATCH --error={error_file}
#SBATCH --time={slurm_args.get('time', '01:00:00')}
#SBATCH --mem-per-cpu={slurm_args.get('memory', '2G')}
#SBATCH --cpus-per-task={slurm_args.get('cpus', 1)}
#SBATCH --partition={slurm_args.get('partition', 'cpu')}"""

    if slurm_args.get('account'):
        script_content += f"\n#SBATCH --account={slurm_args['account']}"
    
    # Add array job support if specified
    if slurm_args.get('array'):
        script_content += f"\n#SBATCH --array={slurm_args['array']}"
    
    script_content += f"""

# Load modules and activate mamba environment
# module load mamba  # Uncomment if needed on your cluster
source ~/.bashrc
mamba activate gensim

# Set working directory
cd {os.getcwd()}

# Print job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Combination: {combo_name}"
echo "Working directory: $(pwd)"
echo ""

# Create output directory
mkdir -p "{combo_dir}"

# Run single combination
echo "Running PLINK simulation for combination {combination_id}:"
echo "  Cohort size: {combination['cohort_size']} (Cases: {combination['num_cases']}, Controls: {combination['num_controls']})"
echo "  Prevalence: {combination['prevalence']:.3f}"
echo "  Total SNPs: {combination['total_snps']} (Causal: {combination['causal_snps']}, Null: {combination['null_snps']})"
echo ""

# Execute Python command
python -c "
import os
from gensim import PLINKSimulationConfig, PLINKSimulationSet, PLINKSimulator

# Create SNP sets
snp_sets = [
    PLINKSimulationSet(
        num_snps={combination['null_snps']},
        label='null',
        min_freq={base_config.min_freq},
        max_freq={base_config.max_freq},
        het_odds_ratio=1.0,
        hom_odds_ratio=1.0
    )
]

# Add causal SNPs if any
if {combination['causal_snps']} > 0:
    snp_sets.append(
        PLINKSimulationSet(
            num_snps={combination['causal_snps']},
            label='causal',
            min_freq={base_config.min_freq},
            max_freq={base_config.max_freq},
            het_odds_ratio={base_config.het_odds_ratio},
            hom_odds_ratio={hom_or_value}
        )
    )

# Create PLINK configuration
config = PLINKSimulationConfig(
    output_prefix='{combo_name}',
    output_dir='{combo_dir}',
    num_cases={combination['num_cases']},
    num_controls={combination['num_controls']},
    disease_prevalence={combination['prevalence']},
    snp_sets=snp_sets,
    plink_executable='{base_config.plink_executable}',
    random_seed={random_seed_value}
)

# Run simulation
simulator = PLINKSimulator(config)
result = simulator.run_simulation()

if result['success']:
    print('✓ Simulation completed successfully!')
    if 'statistics' in result:
        stats = result['statistics']
        print('  Generated: {{}} individuals, {{}} SNPs'.format(stats['total_individuals'], stats['total_snps']))
        print('  Cases: {{}}, Controls: {{}}'.format(stats['cases'], stats['controls']))
else:
    print('✗ Simulation failed!')
    print('  Error: {{}}'.format(result.get('error', 'Unknown error')))
    exit(1)
"

echo ""
echo "Job finished at: $(date)"
"""
    
    # Write script file
    script_file = os.path.join(slurm_output_dir, f"{job_name}.sh")
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_file, 0o755)
    return script_file


def submit_individual_jobs(
    grid_config: PLINKParameterGrid,
    slurm_args: Dict[str, Any],
    max_jobs: int = None,
    dry_run: bool = False
) -> List[str]:
    """
    Submit individual SLURM jobs for each parameter combination.
    
    Args:
        grid_config: Parameter grid configuration
        slurm_args: SLURM job arguments
        max_jobs: Maximum number of jobs to submit (None for all)
        dry_run: If True, create scripts but don't submit jobs
        
    Returns:
        List of job IDs (empty if dry_run=True)
    """
    combinations = grid_config.get_parameter_combinations()
    
    if max_jobs and max_jobs < len(combinations):
        print(f"Limiting to first {max_jobs} combinations (out of {len(combinations)} total)")
        combinations = combinations[:max_jobs]
    
    job_ids = []
    
    print(f"Creating SLURM jobs for {len(combinations)} combinations...")
    
    for i, combination in enumerate(combinations, 1):
        combo_name = grid_config.get_combination_name(combination)
        
        try:
            # Create SLURM script
            script_file = create_single_combination_script(
                combination, i, grid_config, slurm_args
            )
            
            if dry_run:
                print(f"  {i:3d}. Created script: {os.path.basename(script_file)}")
            else:
                # Submit job
                result = subprocess.run(
                    ['sbatch', script_file],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Extract job ID
                job_id = result.stdout.strip().split()[-1]
                job_ids.append(job_id)
                
                print(f"  {i:3d}. Submitted job {job_id}: {combo_name}")
                
        except subprocess.CalledProcessError as e:
            print(f"  {i:3d}. Error submitting {combo_name}: {e}")
            print(f"       STDERR: {e.stderr}")
        except Exception as e:
            print(f"  {i:3d}. Unexpected error for {combo_name}: {e}")
    
    return job_ids


def create_array_job_script(
    grid_config: PLINKParameterGrid,
    slurm_args: Dict[str, Any],
    max_array_size: int = 1000
) -> str:
    """
    Create a SLURM array job script for parameter combinations.
    
    Args:
        grid_config: Parameter grid configuration
        slurm_args: SLURM job arguments
        max_array_size: Maximum array size
        
    Returns:
        Path to created array job script
    """
    combinations = grid_config.get_parameter_combinations()
    total_combinations = len(combinations)
    
    if total_combinations > max_array_size:
        print(f"Warning: {total_combinations} combinations exceeds max array size {max_array_size}")
        print(f"Consider using chunked jobs instead.")
    
    # Create combinations parameter file
    slurm_output_dir = os.path.join(grid_config.grid_output_dir, "slurm_array")
    os.makedirs(slurm_output_dir, exist_ok=True)
    
    combinations_file = os.path.join(slurm_output_dir, "combinations.json")
    with open(combinations_file, 'w') as f:
        json.dump(combinations, f, indent=2)
    
    # Create array job script
    array_script = os.path.join(slurm_output_dir, "array_job.sh")
    
    # Prepare hom_odds_ratio value
    hom_or_value = "'mult'" if grid_config.hom_odds_ratio == 'mult' else str(grid_config.hom_odds_ratio)
    random_seed_value = str(grid_config.random_seed) if grid_config.random_seed else 'None'
    
    # Build the Python script content separately to avoid f-string issues
    python_script = f'''import json
import sys
from gensim import PLINKSimulationConfig, PLINKSimulationSet, PLINKSimulator

# Load combinations
with open('{combinations_file}', 'r') as f:
    combinations = json.load(f)

# Get combination for this array task (1-indexed)
task_id = int('$SLURM_ARRAY_TASK_ID')
if task_id > len(combinations):
    print(f'Array task ID {{task_id}} exceeds number of combinations {{len(combinations)}}')
    sys.exit(1)

combination = combinations[task_id - 1]  # Convert to 0-indexed

# Create combination name  
base_prefix = '{grid_config.base_prefix}'
combo_name = f'{{base_prefix}}_n{{combination["cohort_size"]}}_prev{{combination["prevalence"]:.3f}}_snps{{combination["total_snps"]}}_causal{{combination["causal_snps"]}}'
combo_dir = '{grid_config.grid_output_dir}/' + combo_name

print(f'Processing combination {{task_id}}: {{combo_name}}')
print(f'  Cohort: {{combination["cohort_size"]}} ({{combination["num_cases"]}} cases, {{combination["num_controls"]}} controls)')
print(f'  Prevalence: {{combination["prevalence"]:.3f}}')
print(f'  SNPs: {{combination["total_snps"]}} ({{combination["causal_snps"]}} causal)')

# Create SNP sets
snp_sets = [
    PLINKSimulationSet(
        num_snps=combination['null_snps'],
        label='null',
        min_freq={grid_config.min_freq},
        max_freq={grid_config.max_freq},
        het_odds_ratio=1.0,
        hom_odds_ratio=1.0
    )
]

if combination['causal_snps'] > 0:
    snp_sets.append(
        PLINKSimulationSet(
            num_snps=combination['causal_snps'],
            label='causal',
            min_freq={grid_config.min_freq},
            max_freq={grid_config.max_freq},
            het_odds_ratio={grid_config.het_odds_ratio},
            hom_odds_ratio={hom_or_value}
        )
    )

# Create configuration
config = PLINKSimulationConfig(
    output_prefix=combo_name,
    output_dir=combo_dir,
    num_cases=combination['num_cases'],
    num_controls=combination['num_controls'],
    disease_prevalence=combination['prevalence'],
    snp_sets=snp_sets,
    plink_executable='{grid_config.plink_executable}',
    random_seed={random_seed_value}
)

# Run simulation
simulator = PLINKSimulator(config)
result = simulator.run_simulation()

if result['success']:
    print('✓ Simulation completed successfully!')
else:
    print('✗ Simulation failed!')
    print(f'Error: {{result.get("error", "Unknown error")}}')
    sys.exit(1)'''

    script_content = f"""#!/bin/bash

#SBATCH --job-name=plink_array
#SBATCH --output={slurm_output_dir}/plink_array_%A_%a.out
#SBATCH --error={slurm_output_dir}/plink_array_%A_%a.err
#SBATCH --time={slurm_args.get('time', '01:00:00')}
#SBATCH --mem-per-cpu={slurm_args.get('memory', '2G')}
#SBATCH --cpus-per-task={slurm_args.get('cpus', 1)}
#SBATCH --partition={slurm_args.get('partition', 'cpu')}
#SBATCH --array=1-{min(total_combinations, max_array_size)}"""

    if slurm_args.get('account'):
        script_content += f"\n#SBATCH --account={slurm_args['account']}"

    script_content += f"""

# Load modules and activate mamba environment
# module load mamba  # Uncomment if needed on your cluster
source ~/.bashrc
mamba activate gensim

# Set working directory
cd {os.getcwd()}

# Print job information
echo "Array job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"

# Run combination for this array index
python -c "{python_script}"

echo "Array task finished at: $(date)"
"""
    
    with open(array_script, 'w') as f:
        f.write(script_content)
    
    os.chmod(array_script, 0o755)
    return array_script


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Submit individual parameter combinations as SLURM jobs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit each combination as separate job
  python submit_individual_slurm.py \\
    --cohort-sizes "500,1000" --prevalences "0.01,0.05" \\
    --total-snps "5000,10000" --causal-snps "50,100"
  
  # Use array job for many combinations
  python submit_individual_slurm.py \\
    --cohort-sizes "500,1000,2000" --prevalences "0.01,0.05,0.1" \\
    --total-snps "5000,10000,20000" --causal-snps "50,100,200" \\
    --use-array
  
  # Limit number of jobs
  python submit_individual_slurm.py \\
    --cohort-sizes "500,1000,2000,5000" --prevalences "0.01,0.05,0.1,0.2" \\
    --total-snps "5000,10000,20000,50000" --causal-snps "50,100,200,500" \\
    --max-jobs 50
        """
    )
    
    # Grid parameters (same as main script)
    parser.add_argument("--cohort-sizes", required=True,
                       help="Comma-separated list of cohort sizes")
    parser.add_argument("--prevalences", required=True,
                       help="Comma-separated list of prevalences")
    parser.add_argument("--total-snps", required=True,
                       help="Comma-separated list of total SNP counts")
    parser.add_argument("--causal-snps", required=True,
                       help="Comma-separated list of causal SNP counts")
    parser.add_argument("--grid-output-dir", default="individual_jobs_grid",
                       help="Output directory for parameter grid")
    parser.add_argument("--base-prefix", default="dataset",
                       help="Base prefix for dataset names")
    
    # Configuration
    parser.add_argument("--case-control-ratio", type=float, default=1.0,
                       help="Case-to-control ratio")
    parser.add_argument("--min-freq", type=float, default=0.01,
                       help="Minimum allele frequency")
    parser.add_argument("--max-freq", type=float, default=0.5,
                       help="Maximum allele frequency")
    parser.add_argument("--het-or", type=float, default=1.5,
                       help="Heterozygote odds ratio")
    parser.add_argument("--hom-or", default="mult",
                       help="Homozygote odds ratio")
    parser.add_argument("--plink-executable", default="plink",
                       help="PLINK executable")
    parser.add_argument("--random-seed", type=int,
                       help="Random seed")
    
    # SLURM parameters
    parser.add_argument("--time", default="01:00:00",
                       help="Job time limit")
    parser.add_argument("--memory", default="2G",
                       help="Memory per job")
    parser.add_argument("--cpus", type=int, default=1,
                       help="CPUs per job")
    parser.add_argument("--partition", default="cpu",
                       help="SLURM partition")
    parser.add_argument("--account",
                       help="SLURM account")
    
    # Job control
    parser.add_argument("--max-jobs", type=int,
                       help="Maximum number of jobs to submit")
    parser.add_argument("--use-array", action="store_true",
                       help="Use SLURM array job instead of individual jobs")
    parser.add_argument("--max-array-size", type=int, default=1000,
                       help="Maximum array size for array jobs")
    parser.add_argument("--dry-run", action="store_true",
                       help="Create scripts but don't submit")
    
    args = parser.parse_args()
    
    # Parse parameters
    def parse_int_list(value_str):
        return [int(x.strip()) for x in value_str.split(',')]
    
    def parse_float_list(value_str):
        return [float(x.strip()) for x in value_str.split(',')]
    
    hom_or = args.hom_or if args.hom_or == "mult" else float(args.hom_or)
    
    # Create grid config
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
        print(f"Error creating grid configuration: {e}")
        return 1
    
    total_combinations = len(grid_config.get_parameter_combinations())
    print(f"Total parameter combinations: {total_combinations}")
    
    slurm_args = {
        'time': args.time,
        'memory': args.memory,
        'cpus': args.cpus,
        'partition': args.partition,
        'account': args.account
    }
    
    if args.use_array:
        # Create array job
        print(f"Creating SLURM array job...")
        array_script = create_array_job_script(grid_config, slurm_args, args.max_array_size)
        
        if args.dry_run:
            print(f"Array job script created: {array_script}")
        else:
            try:
                result = subprocess.run(['sbatch', array_script], capture_output=True, text=True, check=True)
                job_id = result.stdout.strip().split()[-1]
                print(f"Submitted array job: {job_id}")
                print(f"Monitor with: squeue -j {job_id}")
            except subprocess.CalledProcessError as e:
                print(f"Error submitting array job: {e}")
                return 1
    else:
        # Submit individual jobs
        job_ids = submit_individual_jobs(grid_config, slurm_args, args.max_jobs, args.dry_run)
        
        if job_ids:
            print(f"\\nSubmitted {len(job_ids)} individual jobs")
            print(f"Monitor with: squeue -u $USER")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
