# Gensim: Genomic Data Simulation with GCTA

A Python package for generating simulated genomic data using GCTA (Genome-wide Complex Trait Analysis) with support for parameter grids and batch processing.

## Features

- **Parameter Grid Support**: Run simulations across multiple combinations of cohort sizes, numbers of causal SNPs, heritabilities, and prevalences
- **Both Trait Types**: Support for quantitative and binary (case-control) traits
- **Flexible Configuration**: JSON-based configuration files or command-line parameters
- **Automated File Management**: Automatic generation of causal SNP lists and individual keep files
- **Comprehensive Logging**: Detailed logging of all simulation steps and results
- **Result Validation**: Automatic validation of output files and generation of summary reports
- **Reproducible Results**: Support for random seeds to ensure reproducibility

## Installation

### Prerequisites

1. **GCTA**: Install GCTA from [https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)
2. **PLINK Binary Files**: You need `.bed`, `.bim`, and `.fam` files for your reference dataset

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

### Command Line Interface

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
- GCTA (external dependency)

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
