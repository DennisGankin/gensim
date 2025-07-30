#!/bin/bash
#SBATCH --job-name=bed_to_h5_conservative
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=bed_to_h5_conservative_%j.out
#SBATCH --error=bed_to_h5_conservative_%j.err

# Load bashrc and environment
source ~/.bashrc

# Print job information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo "Allocated memory: $SLURM_MEM_PER_NODE MB"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"

# Navigate to the script directory
cd /cluster/project/beltrao/gankin/vnn/gensim
mamba activate gensim

# Force chunked processing with very conservative memory usage (8GB chunks)
# This should work even for the largest datasets
python convert_bed_to_h5.py --input-dir="data/ukb_sim" --output-dir="data/ukb_sim_h5py" --chunked --memory-usage 8.0 --verbose

echo "Job completed at: $(date)"
