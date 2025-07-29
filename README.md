# Gensim: Genomic Data Simulation with GCTA and PLINK

A Python package for generating simulated genomic data using GCTA (Genome-wide Complex Trait Analysis) and PLINK with support for parameter grids and batch processing.

## Features

- **Parameter Grid Support**: Run simulations across multiple combinations of cohort sizes, numbers of causal SNPs, heritabilities, and prevalences
- **Both Trait Types**: Support for quantitative and binary (case-control) traits
- **PLINK Dataset Creation**: Generate synthetic SNP datasets using PLINK's simulation capabilities
- **Flexible Configuration**: JSON-based configuration files or command-line parameters
- **Automated File Management**: Automatic generation of causal SNP lists and individual keep files
- **Comprehensive Logging**: Detailed logging of all simulation steps and results
- **Result Validation**: Automatic validation of output files and generation of summary reports
- **Reproducible Results**: Support for random seeds to ensure reproducibility

## Installation

### Prerequisites

1. **GCTA**: Install GCTA from [https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)
2. **PLINK**: Install PLINK from [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/) (for dataset creation)
3. **PLINK Binary Files**: You need `.bed`, `.bim`, and `.fam` files for your reference dataset (for GCTA simulations)

### Install Gensim

```bash
# Clone the repository
git clone https://github.com/DennisGankin/gensim.git
cd gensim

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

## Quick Start

### PLINK Dataset Creation

```bash
# Create a simple synthetic SNP dataset
python main.py --create-plink-dataset --plink-output my_data --plink-cases 1000 --plink-controls 1000

# Create custom dataset with specific parameters
python main.py --create-plink-dataset --plink-output custom_data \
  --plink-cases 2000 --plink-controls 3000 \
  --plink-null-snps 50000 --plink-disease-snps 200 \
  --plink-prevalence 0.05
```

### GCTA Phenotype Simulation

```bash
# Create an example configuration file
python main.py --create-config config.json

# Run simulation grid from config file
python main.py --config config.json

# Run quick simulation with inline parameters
python main.py --bfile test --cohort-sizes 1000 2000 --heritabilities 0.5 0.8

# Run single custom simulation
python main.py --bfile test --single --cohort-size 1000 --num-causal 100 --heritability 0.5
```

### Python API

```python
from gensim import GCTASimulator, SimulationConfig

# Configure simulation parameters
config = SimulationConfig(
    bfile="test",  # Your PLINK binary file prefix
    output_dir="simulations",
    cohort_sizes=[500, 1000, 2000],
    num_causal_snps=[50, 100, 200],
    heritabilities=[0.3, 0.5, 0.8],
    trait_type="quantitative",
    random_seed=42
)

# Run simulations
simulator = GCTASimulator(config)
results = simulator.run_simulation_grid()

# Get summary
summary = simulator.get_results_summary()
print(f"Completed {summary['successful']}/{summary['total_simulations']} simulations")
```

## Configuration

### Example Configuration File

```json
{
  "bfile": "test",
  "output_dir": "simulations",
  "cohort_sizes": [500, 1000, 2000],
  "num_causal_snps": [50, 100, 200],
  "heritabilities": [0.3, 0.5, 0.8],
  "prevalences": [0.05, 0.1, 0.2],
  "num_replications": 1,
  "trait_type": "quantitative",
  "gcta_executable": "gcta64",
  "random_seed": 42
}
```

### Configuration Parameters

- **bfile**: Base name for PLINK binary files (required)
- **output_dir**: Directory for simulation outputs
- **cohort_sizes**: List of sample sizes to simulate
- **num_causal_snps**: List of numbers of causal SNPs
- **heritabilities**: List of heritability values (0-1)
- **prevalences**: List of disease prevalences (for binary traits)
- **trait_type**: "quantitative" or "binary"
- **num_replications**: Number of replications per parameter combination
- **gcta_executable**: Path to GCTA executable
- **random_seed**: Seed for reproducibility
- **causal_snplist**: Pre-defined causal SNP file (optional)
- **keep_individuals**: Pre-defined individuals to keep (optional)

## Examples

### Quantitative Trait Simulation

```python
from gensim import GCTASimulator, SimulationConfig

config = SimulationConfig(
    bfile="test",
    output_dir="quantitative_sims",
    cohort_sizes=[1000, 2000],
    num_causal_snps=[100, 200],
    heritabilities=[0.3, 0.5, 0.8],
    trait_type="quantitative",
    num_replications=3,
    random_seed=42
)

simulator = GCTASimulator(config)
results = simulator.run_simulation_grid()
```

### Binary Trait (Case-Control) Simulation

```python
config = SimulationConfig(
    bfile="test",
    output_dir="binary_sims",
    cohort_sizes=[1000, 2000],
    num_causal_snps=[100, 200],
    heritabilities=[0.5, 0.8],
    prevalences=[0.05, 0.1, 0.2],
    trait_type="binary",
    random_seed=42
)

simulator = GCTASimulator(config)
results = simulator.run_simulation_grid()
```

### Single Custom Simulation

```python
config = SimulationConfig(bfile="test", trait_type="quantitative")
simulator = GCTASimulator(config)

result = simulator.run_single_custom_simulation(
    cohort_size=1000,
    num_causal=100,
    heritability=0.5
)
```

## PLINK Dataset Creation

### Simple Dataset Creation

```python
from gensim import PLINKSimulator, PLINKSimulationConfig, PLINKSimulationSet

# Create simple dataset with null and disease SNPs
snp_sets = PLINKSimulator.create_simple_simulation_sets(
    num_null=10000,     # 10,000 null SNPs
    num_disease=100     # 100 disease-associated SNPs
)

config = PLINKSimulationConfig(
    output_prefix="my_dataset",
    output_dir="plink_data",
    num_cases=1000,
    num_controls=1000,
    disease_prevalence=0.05,
    snp_sets=snp_sets,
    random_seed=42
)

simulator = PLINKSimulator(config)
result = simulator.run_simulation()
```

### Custom SNP Sets

```python
# Create custom SNP sets with different frequency ranges
snp_sets = [
    # Rare variants (MAF < 5%)
    PLINKSimulationSet(5000, "rare_null", 0.001, 0.05, 1.00, 1.00),
    PLINKSimulationSet(20, "rare_disease", 0.001, 0.05, 2.50, "mult"),
    
    # Common variants (MAF 5-50%)
    PLINKSimulationSet(10000, "common_null", 0.05, 0.50, 1.00, 1.00),
    PLINKSimulationSet(80, "common_disease", 0.05, 0.50, 1.80, "mult"),
]

config = PLINKSimulationConfig(
    output_prefix="custom_dataset",
    num_cases=2000,
    num_controls=3000,
    snp_sets=snp_sets
)

simulator = PLINKSimulator(config)
result = simulator.run_simulation()
```

### PLINK Simulation Parameters

Each `PLINKSimulationSet` specifies:
- **num_snps**: Number of SNPs in this set
- **label**: Unique label for SNP naming
- **min_freq**: Minimum allele frequency
- **max_freq**: Maximum allele frequency  
- **het_odds_ratio**: Odds ratio for heterozygotes
- **hom_odds_ratio**: Odds ratio for homozygotes (or "mult" for multiplicative)

### Command Line PLINK Usage

```bash
# Create simple dataset
python main.py --create-plink-dataset --plink-output my_data \
  --plink-cases 1000 --plink-controls 1000

# Create custom dataset
python main.py --create-plink-dataset --plink-output gwas_data \
  --plink-cases 5000 --plink-controls 5000 \
  --plink-null-snps 100000 --plink-disease-snps 500 \
  --plink-prevalence 0.01
```

## Parameter Grid Simulation

### Local Parameter Grid

Create multiple datasets with different parameter combinations:

```bash
# Create parameter grid locally
python main.py --create-parameter-grid \
  --grid-cohort-sizes "500,1000,2000" \
  --grid-prevalences "0.01,0.05,0.1" \
  --grid-total-snps "5000,10000,20000" \
  --grid-causal-snps "50,100,200" \
  --grid-output-dir "my_parameter_grid"
```

This creates 3×3×3×3 = 81 different dataset combinations.

### SLURM Cluster Parallel Execution

For large parameter grids, submit jobs to SLURM clusters for parallel processing:

#### Chunked Jobs (Recommended for medium grids)

```bash
# Submit parameter grid as chunked SLURM jobs
python submit_parameter_grid_slurm.py \
  --cohort-sizes "500,1000,2000,5000" \
  --prevalences "0.001,0.01,0.05,0.1,0.2" \
  --total-snps "5000,10000,20000,50000" \
  --causal-snps "50,100,200,500" \
  --num-jobs 16 \
  --time "04:00:00" \
  --memory "8G" \
  --partition "cpu"
```

This splits 4×5×4×4 = 320 combinations into 16 SLURM jobs (~20 combinations each).

#### Individual Jobs (Maximum parallelization)

```bash
# Submit each combination as separate SLURM job
python submit_individual_slurm.py \
  --cohort-sizes "1000,2000,5000" \
  --prevalences "0.01,0.05,0.1" \
  --total-snps "10000,50000" \
  --causal-snps "100,500" \
  --time "02:00:00" \
  --memory "4G" \
  --partition "cpu"
```

Creates 3×3×2×2 = 36 individual SLURM jobs.

#### Array Jobs (Most efficient for large grids)

```bash
# Use SLURM array jobs for efficient submission
python submit_individual_slurm.py \
  --cohort-sizes "500,1000,2000,5000,10000" \
  --prevalences "0.001,0.01,0.05,0.1,0.2,0.5" \
  --total-snps "5000,10000,20000,50000,100000" \
  --causal-snps "50,100,200,500,1000" \
  --use-array \
  --time "03:00:00" \
  --memory "6G" \
  --partition "cpu"
```

Creates 1 array job with 5×6×5×5 = 750 tasks.

#### SLURM Job Monitoring

Monitor submitted jobs:

```bash
# Check job status
squeue -u $USER

# Monitor specific job
squeue -j 12345

# Use automatic monitoring script (created by chunked submission)
./parameter_grid_output/monitor_jobs.sh

# Check results
find parameter_grid_output -name "grid_summary.txt"
```

#### Custom SLURM Parameters

Customize resource requirements for your cluster:

```bash
python submit_parameter_grid_slurm.py \
  --cohort-sizes "10000,50000" \
  --prevalences "0.01,0.1" \
  --total-snps "500000,1000000" \
  --causal-snps "1000,5000" \
  --num-jobs 8 \
  --time "12:00:00" \        # 12 hours
  --memory "32G" \           # 32GB RAM
  --cpus "4" \               # 4 CPUs per job
  --partition "highmem" \    # High-memory partition
  --account "myproject"      # Billing account
```

#### SLURM Examples and Testing

```bash
# View comprehensive examples
python slurm_examples.py

# Test SLURM submission (dry run)
python test_slurm_submission.py

# Test with actual job submission
python test_slurm_submission.py --submit
```

## Output Files

For each simulation, the following files are generated:

- **{name}.par**: Parameter file with simulation details
- **{name}.phen**: Phenotype file with simulated trait values
- **{name}.log**: GCTA log file
- **simulation_summary.csv**: Summary of all simulations in the batch

## GCTA Command Example

The package generates commands like:

```bash
gcta64 --bfile test --simu-qt --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-rep 3 --keep test.indi.list --out test
```

## Dependencies

- Python ≥ 3.7
- pandas ≥ 1.0.0
- numpy ≥ 1.18.0
- GCTA (external dependency for phenotype simulation)
- PLINK (external dependency for dataset creation)

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
Gankin, D. (2025). Gensim: A Python package for genomic data simulation using GCTA. 
GitHub repository: https://github.com/DennisGankin/gensim
```

## Support

For questions, issues, or feature requests, please open an issue on GitHub.
