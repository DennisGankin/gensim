#!/bin/bash

# Script to run GCTA simulations for multiple parameter combinations
# Usage: ./run_gcta_simulations.sh

# Define parameter arrays
COHORT_SIZES=(5000 50000 500000)
PREVALENCES=(0.01 0.05 0.10 0.20 0.50)
HERITABILITIES=(0.1 0.3 0.5 0.7)
CAUSAL_SNPS=(80 320 1280 5120)

# Fixed parameters
GENOTYPE_BASE_DIR="data/ukb_plink"
CAUSAL_SNPS_DIR="data/ukb_plink"
OUTPUT_BASE="data/ukb_plink/gcta"
SLURM_SCRIPT_DIR="slurm_scripts"

# Create necessary directories
mkdir -p "$OUTPUT_BASE"
mkdir -p "$SLURM_SCRIPT_DIR"

# Function to calculate cases and controls from cohort size and prevalence
calculate_cases_controls() {
    local cohort_size=$1
    local prevalence=$2
    
    # Calculate number of cases (round to nearest integer)
    local cases=$(echo "$cohort_size * $prevalence" | bc -l | xargs printf "%.0f")
    
    # Ensure at least 1 case and at least 1 control
    if [ "$cases" -eq 0 ]; then
        cases=1
    fi
    if [ "$cases" -eq "$cohort_size" ]; then
        cases=$((cohort_size - 1))
    fi
    
    local controls=$((cohort_size - cases))
    
    echo "$cases $controls"
}

# Function to create SLURM script for a single simulation
create_slurm_script() {
    local cohort_size=$1
    local prevalence=$2
    local heritability=$3
    local causal_snps=$4
    local cases=$5
    local controls=$6
    
    # Create descriptive job name and output prefix
    local job_name="gcta_n${cohort_size}_prev${prevalence}_h2${heritability}_causal${causal_snps}"
    local output_prefix="${OUTPUT_BASE}/${job_name}"
    local slurm_file="${SLURM_SCRIPT_DIR}/${job_name}.slurm"
    local causal_file="${CAUSAL_SNPS_DIR}/causal_snps_${causal_snps}_${cohort_size}.txt"
    local genotype_prefix="${GENOTYPE_BASE_DIR}/genotype_${cohort_size}"
    
    # Create SLURM script
    cat > "$slurm_file" << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_scripts/${job_name}_%j.out
#SBATCH --error=slurm_scripts/${job_name}_%j.err

# Load environment
source ~/.bashrc
mamba activate gensim

# Navigate to the project directory

# Print job information
echo "Job started at: \$(date)"
echo "Running on node: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Parameters:"
echo "  Cohort size: ${cohort_size}"
echo "  Prevalence: ${prevalence}"
echo "  Heritability: ${heritability}"
echo "  Causal SNPs: ${causal_snps}"
echo "  Cases: ${cases}"
echo "  Controls: ${controls}"
echo "  Causal SNPs file: ${causal_file}"
echo "  Genotype prefix: ${genotype_prefix}"
echo "  Output prefix: ${output_prefix}"

# Check if causal SNPs file exists
if [ ! -f "${causal_file}" ]; then
    echo "ERROR: Causal SNPs file not found: ${causal_file}"
    echo "Please run select_causal_snps.sh first to generate the causal SNPs files"
    exit 1
fi

# Check if genotype file exists
if [ ! -f "${genotype_prefix}.bed" ]; then
    echo "ERROR: Genotype file not found: ${genotype_prefix}.bed"
    echo "Please ensure the genotype files are available"
    exit 1
fi

# Create output directory
mkdir -p "\$(dirname "${output_prefix}")"

# Run GCTA simulation
echo "Running GCTA simulation..."
gcta \\
    --bfile ${genotype_prefix} \\
    --simu-cc ${cases} ${controls} \\
    --simu-hsq ${heritability} \\
    --simu-causal-loci ${causal_file} \\
    --out ${output_prefix}

# Check if simulation was successful
if [ \$? -eq 0 ]; then
    echo "GCTA simulation completed successfully"
    echo "Output files:"
    ls -la ${output_prefix}.*
else
    echo "ERROR: GCTA simulation failed"
    exit 1
fi

echo "Job completed at: \$(date)"
EOF

    echo "$slurm_file"
}

# Function to submit or save all jobs
submit_jobs() {
    local submit_mode=$1  # "submit" or "save"
    local job_count=0
    local submitted_jobs=()
    
    echo "Processing parameter combinations..."
    echo "Mode: $submit_mode"
    
    for cohort_size in "${COHORT_SIZES[@]}"; do
        for prevalence in "${PREVALENCES[@]}"; do
            for heritability in "${HERITABILITIES[@]}"; do
                for causal_snps in "${CAUSAL_SNPS[@]}"; do
                    # Calculate cases and controls
                    cases_controls=$(calculate_cases_controls "$cohort_size" "$prevalence")
                    cases=$(echo "$cases_controls" | cut -d' ' -f1)
                    controls=$(echo "$cases_controls" | cut -d' ' -f2)
                    
                    # Create SLURM script
                    slurm_file=$(create_slurm_script "$cohort_size" "$prevalence" "$heritability" "$causal_snps" "$cases" "$controls")
                    
                    job_count=$((job_count + 1))
                    
                    if [ "$submit_mode" = "submit" ]; then
                        # Submit the job
                        job_id=$(sbatch "$slurm_file" | grep -o '[0-9]*')
                        submitted_jobs+=("$job_id")
                        echo "[$job_count] Submitted job $job_id: $(basename "$slurm_file")"
                    else
                        echo "[$job_count] Created: $slurm_file"
                    fi
                done
            done
        done
    done
    
    echo ""
    echo "Summary:"
    echo "  Total parameter combinations: $job_count"
    
    if [ "$submit_mode" = "submit" ]; then
        echo "  Submitted jobs: ${#submitted_jobs[@]}"
        echo "  Job IDs: ${submitted_jobs[*]}"
        echo ""
        echo "Monitor jobs with: squeue -u \$USER"
        echo "Cancel all jobs with: scancel ${submitted_jobs[*]}"
    else
        echo "  SLURM scripts created in: $SLURM_SCRIPT_DIR"
        echo ""
        echo "To submit all jobs, run:"
        echo "  for script in $SLURM_SCRIPT_DIR/*.slurm; do sbatch \"\$script\"; done"
    fi
}

# Main execution
echo "GCTA Simulation Parameter Grid"
echo "=============================="
echo ""
echo "Parameter ranges:"
echo "  Cohort sizes: ${COHORT_SIZES[*]}"
echo "  Prevalences: ${PREVALENCES[*]}"
echo "  Heritabilities: ${HERITABILITIES[*]}"
echo "  Causal SNPs: ${CAUSAL_SNPS[*]}"
echo ""

# Calculate total combinations
total_combinations=$((${#COHORT_SIZES[@]} * ${#PREVALENCES[@]} * ${#HERITABILITIES[@]} * ${#CAUSAL_SNPS[@]}))
echo "Total parameter combinations: $total_combinations"
echo ""

# Check if user wants to just create scripts or submit jobs (default is submit)
if [ "$1" = "--dry-run" ]; then
    echo "Creating SLURM scripts only (use without --dry-run to submit jobs)..."
    submit_jobs "save"
else
    echo "Will submit all jobs to SLURM..."
    read -p "Continue? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        submit_jobs "submit"
    else
        echo "Cancelled."
        exit 0
    fi
fi

echo ""
echo "Done!"
