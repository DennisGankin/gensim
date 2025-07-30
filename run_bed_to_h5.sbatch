#!/bin/bash
#SBATCH --job-name=bed_to_h5_conversion
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=32
#SBATCH --output=data/bed_to_h5_%j.out
#SBATCH --error=data/bed_to_h5_%j.err

# Load bashrc and environment
source ~/.bashrc

# Activate conda environment if needed (uncomment and modify as needed)
# conda activate your_environment_name

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

# Run the conversion with auto mode (recommended)
python convert_bed_to_h5.py --input-dir="data/ukb_sim" --output-dir="data/ukb_sim_h5py" --auto --verbose

echo "Job completed at: $(date)"
