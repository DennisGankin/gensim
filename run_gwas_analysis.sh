#!/bin/bash
# Script to run PLINK GWAS analysis for multiple parameter combinations
# Usage: ./run_gwas_analysis.sh

# Define parameter arrays (should match the GCTA simulation parameters)# Calculate total combinations
#total_combinations=$((${#COHORT_SIZES[@]} * ${#PREVALENCES[@]} * ${#HERITABILITIES[@]} * ${#CAUSAL_SNPS[@]}))OHORT_SIZES=(500000)
COHORT_SIZES=(5000 50000 500000)
PREVALENCES=(0.01 0.05 0.10 0.20 0.50)
HERITABILITIES=(0.1 0.3 0.5 0.7)
CAUSAL_SNPS=(80 320 1280 5120)

# Fixed parameters
GENOTYPE_BASE_DIR="data/ukb_plink"
GCTA_OUTPUT_DIR="data/ukb_plink/gcta"
GWAS_OUTPUT_DIR="data/results/gwas_results"
SLURM_SCRIPT_DIR="slurm_scripts/gwas"

# Create necessary directories
mkdir -p "$GWAS_OUTPUT_DIR"
mkdir -p "$SLURM_SCRIPT_DIR"

# Function to create SLURM script for a single GWAS analysis
create_gwas_slurm_script() {
    local cohort_size=$1
    local prevalence=$2
    local heritability=$3
    local causal_snps=$4
    
    # Create descriptive job name and file paths
    local job_name="gwas_n${cohort_size}_prev${prevalence}_h2${heritability}_causal${causal_snps}"
    local slurm_file="${SLURM_SCRIPT_DIR}/${job_name}.slurm"
    local genotype_prefix="${GENOTYPE_BASE_DIR}/genotype_${cohort_size}"
    local gcta_prefix="gcta_n${cohort_size}_prev${prevalence}_h2${heritability}_causal${causal_snps}"
    local pheno_file="${GCTA_OUTPUT_DIR}/${gcta_prefix}.phen"
    local gwas_output_prefix="${GWAS_OUTPUT_DIR}/${job_name}"
    
    # Create SLURM script
    cat > "$slurm_file" << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_scripts_gwas/${job_name}_%j.out
#SBATCH --error=slurm_scripts_gwas/${job_name}_%j.err

# Load environment
source ~/.bashrc
mamba activate gensim

# Navigate to the project directory


# Print job information
echo "Job started at: \$(date)"
echo "Running on node: \$(hostname)"
echo "Job ID: \$SLURM_JOB_ID"
echo "GWAS Analysis Parameters:"
echo "  Cohort size: ${cohort_size}"
echo "  Prevalence: ${prevalence}"
echo "  Heritability: ${heritability}"
echo "  Causal SNPs: ${causal_snps}"
echo "  Genotype prefix: ${genotype_prefix}"
echo "  Phenotype file: ${pheno_file}"
echo "  Output prefix: ${gwas_output_prefix}"

# Check if genotype file exists
if [ ! -f "${genotype_prefix}.bed" ]; then
    echo "ERROR: Genotype file not found: ${genotype_prefix}.bed"
    echo "Please ensure the genotype files are available"
    exit 1
fi

# Check if phenotype file exists (from GCTA simulation)
if [ ! -f "${pheno_file}" ]; then
    echo "ERROR: Phenotype file not found: ${pheno_file}"
    echo "Please run GCTA simulations first to generate phenotype files"
    exit 1
fi

# Create output directory
mkdir -p "\$(dirname "${gwas_output_prefix}")"

# Run PLINK GWAS analysis
echo "Running PLINK GWAS analysis..."
plink \\
    --bfile ${genotype_prefix} \\
    --pheno ${pheno_file} \\
    --assoc \\
    --allow-no-sex \\
    --out ${gwas_output_prefix}

# Check if GWAS was successful
if [ \$? -eq 0 ]; then
    echo "PLINK GWAS analysis completed successfully"
    echo "Output files:"
    ls -la ${gwas_output_prefix}.*
    
    # Show summary statistics
    if [ -f "${gwas_output_prefix}.assoc" ]; then
        echo ""
        echo "GWAS Results Summary:"
        echo "Total SNPs analyzed: \$(wc -l < ${gwas_output_prefix}.assoc)"
        echo "Top 10 most significant associations:"
        head -1 ${gwas_output_prefix}.assoc
        tail -n +2 ${gwas_output_prefix}.assoc | sort -k9,9g | head -10
    fi
else
    echo "ERROR: PLINK GWAS analysis failed"
    exit 1
fi

echo "Job completed at: \$(date)"
EOF

    echo "$slurm_file"
}

# Function to submit or save all GWAS jobs
submit_gwas_jobs() {
    local submit_mode=$1  # "submit" or "save"
    local job_count=0
    local submitted_jobs=()
    local skipped_jobs=0
    
    echo "Processing GWAS parameter combinations..."
    echo "Mode: $submit_mode"
    
    for cohort_size in "${COHORT_SIZES[@]}"; do
        for prevalence in "${PREVALENCES[@]}"; do
            for heritability in "${HERITABILITIES[@]}"; do
                for causal_snps in "${CAUSAL_SNPS[@]}"; do
                    # Check if corresponding GCTA output exists
                    local gcta_prefix="gcta_n${cohort_size}_prev${prevalence}_h2${heritability}_causal${causal_snps}"
                    local pheno_file="${GCTA_OUTPUT_DIR}/${gcta_prefix}.phen"
                    
                    if [ ! -f "$pheno_file" ]; then
                        echo "Skipping: $gcta_prefix (phenotype file not found)"
                        skipped_jobs=$((skipped_jobs + 1))
                        continue
                    fi
                    
                    # Create SLURM script
                    slurm_file=$(create_gwas_slurm_script "$cohort_size" "$prevalence" "$heritability" "$causal_snps")
                    
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
    echo "  Total parameter combinations: $((job_count + skipped_jobs))"
    echo "  GWAS jobs processed: $job_count"
    echo "  Skipped (missing phenotype files): $skipped_jobs"
    
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
    
    if [ $skipped_jobs -gt 0 ]; then
        echo ""
        echo "Note: Run GCTA simulations first to generate missing phenotype files"
    fi
}

# Function to check GCTA simulation status
check_gcta_status() {
    local total_expected=0
    local found_files=0
    
    echo "Checking GCTA simulation status..."
    echo "Expected phenotype files:"
    
    for cohort_size in "${COHORT_SIZES[@]}"; do
        for prevalence in "${PREVALENCES[@]}"; do
            for heritability in "${HERITABILITIES[@]}"; do
                for causal_snps in "${CAUSAL_SNPS[@]}"; do
                    local gcta_prefix="gcta_n${cohort_size}_prev${prevalence}_h2${heritability}_causal${causal_snps}"
                    local pheno_file="${GCTA_OUTPUT_DIR}/${gcta_prefix}.phen"
                    
                    total_expected=$((total_expected + 1))
                    
                    if [ -f "$pheno_file" ]; then
                        echo "  ✓ $pheno_file"
                        found_files=$((found_files + 1))
                    else
                        echo "  ✗ $pheno_file"
                    fi
                done
            done
        done
    done
    
    echo ""
    echo "GCTA Status: $found_files/$total_expected phenotype files found"
    
    if [ $found_files -eq 0 ]; then
        echo ""
        echo "No GCTA simulation results found!"
        echo "Please run the GCTA simulation script first:"
        echo "  ./run_gcta_simulations.sh"
        return 1
    elif [ $found_files -lt $total_expected ]; then
        echo ""
        echo "Some GCTA simulations are missing or incomplete."
        echo "GWAS analysis will only run for available phenotype files."
    else
        echo ""
        echo "All GCTA simulations found! Ready for GWAS analysis."
    fi
    
    return 0
}

# Main execution
echo "PLINK GWAS Analysis Parameter Grid"
echo "=================================="
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

# Check GCTA simulation status first
if [ "$1" != "--force" ]; then
    check_gcta_status
    if [ $? -ne 0 ] && [ "$1" != "--dry-run" ]; then
        echo "Use --force to proceed anyway, or --dry-run to create scripts only"
        exit 1
    fi
fi

# Check if user wants to just create scripts or submit jobs (default is submit)
if [ "$1" = "--dry-run" ]; then
    echo "Creating SLURM scripts only (use without --dry-run to submit jobs)..."
    submit_gwas_jobs "save"
elif [ "$1" = "--force" ]; then
    echo "Force mode: Will submit all jobs regardless of missing files..."
    read -p "Continue? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        submit_gwas_jobs "submit"
    else
        echo "Cancelled."
        exit 0
    fi
else
    echo "Will submit GWAS jobs to SLURM..."
    read -p "Continue? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        submit_gwas_jobs "submit"
    else
        echo "Cancelled."
        exit 0
    fi
fi

echo ""
echo "Done!"
