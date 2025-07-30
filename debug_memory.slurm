#!/bin/bash
#SBATCH --job-name=debug_memory
#SBATCH --time=10:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=debug_memory_%j.out
#SBATCH --error=debug_memory_%j.err

# Load bashrc and environment
source ~/.bashrc

# Navigate to the script directory
cd /cluster/project/beltrao/gankin/vnn/gensim
mamba activate gensim

# Run debug script
python debug_memory.py
