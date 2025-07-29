#!/usr/bin/env python3
"""
Command-line interface for gensim - genomic data simulation using GCTA.
"""

import argparse
import json
import sys
from pathlib import Path

from gensim import (
    GCTASimulator, 
    SimulationConfig, 
    PLINKSimulator, 
    PLINKSimulationConfig, 
    PLINKSimulationSet,
    PLINKParameterGrid,
    PLINKParameterGridSimulator
)


def create_example_config(output_file: str):
    """Create an example configuration file."""
    example_config = {
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
    
    with open(output_file, 'w') as f:
        json.dump(example_config, f, indent=2)
    
    print(f"Example configuration saved to: {output_file}")


def load_config_from_file(config_file: str) -> SimulationConfig:
    """Load configuration from JSON file."""
    with open(config_file, 'r') as f:
        config_dict = json.load(f)
    
    return SimulationConfig(**config_dict)


def create_plink_dataset(args) -> int:
    """Create synthetic SNP dataset using PLINK."""
    try:
        print("Creating synthetic SNP dataset with PLINK...")
        
        # Create SNP sets
        snp_sets = PLINKSimulator.create_simple_simulation_sets(
            num_null=args.plink_null_snps,
            num_disease=args.plink_disease_snps
        )
        
        # Create PLINK configuration
        plink_config = PLINKSimulationConfig(
            output_prefix=args.plink_output,
            output_dir=args.output_dir,
            num_cases=args.plink_cases,
            num_controls=args.plink_controls,
            disease_prevalence=args.plink_prevalence,
            snp_sets=snp_sets,
            plink_executable=args.plink_executable,
            random_seed=args.random_seed
        )
        
        # Run PLINK simulation
        simulator = PLINKSimulator(plink_config)
        result = simulator.run_simulation()
        
        if result['success']:
            print(f"✓ PLINK dataset created successfully!")
            print(f"  Output prefix: {result['output_prefix']}")
            print(f"  Output files: {result['output_files']}")
            
            if 'statistics' in result:
                stats = result['statistics']
                print(f"  Total SNPs: {stats['total_snps']}")
                print(f"  Total individuals: {stats['total_individuals']}")
                print(f"  Cases: {stats['cases']}")
                print(f"  Controls: {stats['controls']}")
            
            print(f"  Simulation file: {result['simulation_file']}")
            return 0
        else:
            print(f"✗ PLINK dataset creation failed: {result['error']}")
            return 1
            
    except Exception as e:
        print(f"Error creating PLINK dataset: {e}")
        return 1


def create_parameter_grid(args) -> int:
    """Create multiple PLINK datasets with parameter grid combinations."""
    try:
        print("Creating parameter grid of PLINK datasets...")
        
        # Parse grid parameters (support comma-separated lists)
        def parse_int_list(value_str):
            if not value_str:
                return []
            return [int(x.strip()) for x in value_str.split(',')]
        
        def parse_float_list(value_str):
            if not value_str:
                return []
            return [float(x.strip()) for x in value_str.split(',')]
        
        # Parse homozygote odds ratio
        hom_or = args.grid_hom_or
        if hom_or != "mult":
            try:
                hom_or = float(hom_or)
            except ValueError:
                raise ValueError(f"Invalid homozygote odds ratio: {hom_or}. Must be 'mult' or a number.")
        
        # Build parameter grid configuration
        grid_config = PLINKParameterGrid(
            cohort_sizes=parse_int_list(args.grid_cohort_sizes) if args.grid_cohort_sizes else [1000],
            prevalences=parse_float_list(args.grid_prevalences) if args.grid_prevalences else [0.01],
            total_snps=parse_int_list(args.grid_total_snps) if args.grid_total_snps else [10000],
            causal_snps=parse_int_list(args.grid_causal_snps) if args.grid_causal_snps else [100],
            grid_output_dir=args.grid_output_dir,
            base_prefix=args.grid_prefix,
            case_control_ratio=args.grid_case_control_ratio,
            min_freq=args.grid_min_freq,
            max_freq=args.grid_max_freq,
            het_odds_ratio=args.grid_het_or,
            hom_odds_ratio=hom_or,
            plink_executable=args.plink_executable,
            random_seed=args.random_seed
        )
        
        # Run parameter grid simulation
        grid_simulator = PLINKParameterGridSimulator(grid_config)
        results = grid_simulator.run_parameter_grid()
        
        # Print summary
        print(f"\n{'='*60}")
        print(f"PARAMETER GRID SIMULATION COMPLETE")
        print(f"{'='*60}")
        print(f"Total combinations: {results['total_combinations']}")
        print(f"Successful: {results['successful']}")
        print(f"Failed: {results['failed']}")
        print(f"Success rate: {results['successful']/results['total_combinations']*100:.1f}%")
        print(f"Output directory: {grid_config.grid_output_dir}")
        
        if results['failed'] > 0:
            print(f"\nNote: {results['failed']} combinations failed. Check the summary report for details.")
            return 1
        else:
            print(f"\n✓ All datasets created successfully!")
            return 0
            
    except Exception as e:
        print(f"Error creating parameter grid: {e}")
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Generate simulated genomic data using GCTA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create example config file
  python main.py --create-config config.json
  
  # Run simulation grid from config file
  python main.py --config config.json
  
  # Run quick simulation with inline parameters
  python main.py --bfile test --cohort-sizes 1000 2000 --heritabilities 0.5 0.8
  
  # Run single custom simulation
  python main.py --bfile test --single --cohort-size 1000 --num-causal 100 --heritability 0.5
  
  # Create PLINK dataset (synthetic SNP data)
  python main.py --create-plink-dataset --plink-output test_data --plink-cases 1000 --plink-controls 1000
  
  # Create custom PLINK dataset with specific SNP sets
  python main.py --create-plink-dataset --plink-output custom_data --plink-null-snps 50000 --plink-disease-snps 200
  
  # Create parameter grid of PLINK datasets
  python main.py --create-parameter-grid --grid-cohort-sizes "500,1000,2000" --grid-prevalences "0.01,0.05,0.1"
  
  # Create comprehensive parameter grid
  python main.py --create-parameter-grid \\
    --grid-cohort-sizes "500,1000,2000" \\
    --grid-prevalences "0.01,0.05,0.1" \\
    --grid-total-snps "5000,10000,20000" \\
    --grid-causal-snps "50,100,200" \\
    --grid-output-dir "my_parameter_grid"
        """
    )
    
    # Configuration options
    parser.add_argument("--config", 
                       help="JSON configuration file")
    parser.add_argument("--create-config", 
                       help="Create example configuration file and exit")
    
    # Basic required parameters
    parser.add_argument("--bfile", 
                       help="Base name for PLINK binary files")
    parser.add_argument("--output-dir", 
                       default="simulations",
                       help="Output directory for simulations")
    
    # Grid parameters
    parser.add_argument("--cohort-sizes", 
                       type=int, 
                       nargs='+',
                       help="List of cohort sizes")
    parser.add_argument("--num-causal-snps", 
                       type=int, 
                       nargs='+',
                       help="List of numbers of causal SNPs")
    parser.add_argument("--heritabilities", 
                       type=float, 
                       nargs='+',
                       help="List of heritability values")
    parser.add_argument("--prevalences", 
                       type=float, 
                       nargs='+',
                       help="List of prevalence values (for binary traits)")
    
    # Single simulation mode
    parser.add_argument("--single", 
                       action="store_true",
                       help="Run single simulation instead of grid")
    parser.add_argument("--cohort-size", 
                       type=int,
                       help="Cohort size for single simulation")
    parser.add_argument("--num-causal", 
                       type=int,
                       help="Number of causal SNPs for single simulation")
    parser.add_argument("--heritability", 
                       type=float,
                       help="Heritability for single simulation")
    parser.add_argument("--prevalence", 
                       type=float,
                       help="Prevalence for single simulation (binary traits)")
    
    # Other options
    parser.add_argument("--trait-type", 
                       choices=["quantitative", "binary"],
                       default="quantitative",
                       help="Type of trait to simulate")
    parser.add_argument("--num-replications", 
                       type=int, 
                       default=1,
                       help="Number of replications per simulation")
    parser.add_argument("--gcta-executable", 
                       default="gcta64",
                       help="GCTA executable name or path")
    parser.add_argument("--random-seed", 
                       type=int,
                       help="Random seed for reproducibility")
    parser.add_argument("--causal-snplist", 
                       help="File with pre-defined causal SNPs")
    parser.add_argument("--keep-individuals", 
                       help="File with individuals to keep")
    parser.add_argument("--no-cleanup", 
                       action="store_true",
                       help="Don't clean up temporary files")
    
    # PLINK dataset creation options
    parser.add_argument("--create-plink-dataset", 
                       action="store_true",
                       help="Create synthetic SNP dataset using PLINK")
    parser.add_argument("--plink-output", 
                       default="plink_data",
                       help="Output prefix for PLINK dataset")
    parser.add_argument("--plink-cases", 
                       type=int, 
                       default=100,
                       help="Number of cases in PLINK dataset")
    parser.add_argument("--plink-controls", 
                       type=int, 
                       default=100,
                       help="Number of controls in PLINK dataset")
    parser.add_argument("--plink-prevalence", 
                       type=float, 
                       default=0.01,
                       help="Disease prevalence for PLINK simulation")
    parser.add_argument("--plink-null-snps", 
                       type=int, 
                       default=10000,
                       help="Number of null (non-associated) SNPs")
    parser.add_argument("--plink-disease-snps", 
                       type=int, 
                       default=100,
                       help="Number of disease-associated SNPs")
    parser.add_argument("--plink-executable", 
                       default="plink",
                       help="PLINK executable name or path")
    
    # Parameter grid options
    parser.add_argument("--create-parameter-grid", 
                       action="store_true",
                       help="Create parameter grid of PLINK datasets")
    parser.add_argument("--grid-cohort-sizes", 
                       type=str,
                       help="Comma-separated list of cohort sizes (e.g., '500,1000,2000')")
    parser.add_argument("--grid-prevalences", 
                       type=str,
                       help="Comma-separated list of prevalences (e.g., '0.01,0.05,0.1')")
    parser.add_argument("--grid-total-snps", 
                       type=str,
                       help="Comma-separated list of total SNP counts (e.g., '5000,10000,20000')")
    parser.add_argument("--grid-causal-snps", 
                       type=str,
                       help="Comma-separated list of causal SNP counts (e.g., '50,100,200')")
    parser.add_argument("--grid-output-dir", 
                       default="parameter_grid",
                       help="Output directory for parameter grid")
    parser.add_argument("--grid-prefix", 
                       default="dataset",
                       help="Base prefix for dataset names")
    parser.add_argument("--grid-case-control-ratio", 
                       type=float, 
                       default=1.0,
                       help="Case-to-control ratio (default: 1.0 for 1:1 ratio)")
    parser.add_argument("--grid-min-freq", 
                       type=float, 
                       default=0.01,
                       help="Minimum allele frequency for grid datasets")
    parser.add_argument("--grid-max-freq", 
                       type=float, 
                       default=0.5,
                       help="Maximum allele frequency for grid datasets")
    parser.add_argument("--grid-het-or", 
                       type=float, 
                       default=1.5,
                       help="Heterozygote odds ratio for causal SNPs")
    parser.add_argument("--grid-hom-or", 
                       default="mult",
                       help="Homozygote odds ratio for causal SNPs ('mult' or float)")
    
    args = parser.parse_args()
    
    # Handle create-config option
    if args.create_config:
        create_example_config(args.create_config)
        return 0
    
    # Handle PLINK dataset creation
    if args.create_plink_dataset:
        return create_plink_dataset(args)
    
    # Handle parameter grid creation
    if args.create_parameter_grid:
        return create_parameter_grid(args)
    
    # Load configuration
    if args.config:
        try:
            config = load_config_from_file(args.config)
        except Exception as e:
            print(f"Error loading configuration file: {e}")
            return 1
    else:
        # Build configuration from command line arguments
        if not args.bfile:
            print("Error: --bfile is required when not using --config, --create-plink-dataset, or --create-parameter-grid")
            return 1
        
        config_dict = {
            "bfile": args.bfile,
            "output_dir": args.output_dir,
            "trait_type": args.trait_type,
            "num_replications": args.num_replications,
            "gcta_executable": args.gcta_executable
        }
        
        # Add optional parameters
        if args.random_seed is not None:
            config_dict["random_seed"] = args.random_seed
        if args.causal_snplist:
            config_dict["causal_snplist"] = args.causal_snplist
        if args.keep_individuals:
            config_dict["keep_individuals"] = args.keep_individuals
        
        # Add grid parameters or single simulation parameters
        if args.single:
            # Single simulation mode
            if not all([args.cohort_size, args.num_causal, args.heritability]):
                print("Error: --cohort-size, --num-causal, and --heritability are required for single simulation")
                return 1
            
            config_dict.update({
                "cohort_sizes": [args.cohort_size],
                "num_causal_snps": [args.num_causal],
                "heritabilities": [args.heritability]
            })
            
            if args.prevalence is not None:
                config_dict["prevalences"] = [args.prevalence]
        else:
            # Grid simulation mode
            config_dict.update({
                "cohort_sizes": args.cohort_sizes or [1000],
                "num_causal_snps": args.num_causal_snps or [100],
                "heritabilities": args.heritabilities or [0.5]
            })
            
            if args.prevalences:
                config_dict["prevalences"] = args.prevalences
        
        try:
            config = SimulationConfig(**config_dict)
        except Exception as e:
            print(f"Error creating configuration: {e}")
            return 1
    
    # Run simulations
    try:
        simulator = GCTASimulator(config)
        
        if args.single:
            # Run single simulation
            result = simulator.run_single_custom_simulation(
                cohort_size=config.cohort_sizes[0],
                num_causal=config.num_causal_snps[0],
                heritability=config.heritabilities[0],
                prevalence=config.prevalences[0] if config.prevalences else None
            )
            
            if result['success']:
                print(f"Single simulation completed successfully!")
                print(f"Output files: {result['output_files']}")
            else:
                print(f"Simulation failed: {result['error']}")
                return 1
        else:
            # Run simulation grid
            results = simulator.run_simulation_grid(cleanup_temp=not args.no_cleanup)
            
            # Print summary
            summary = simulator.get_results_summary()
            print(f"\nSimulation Summary:")
            print(f"Total simulations: {summary['total_simulations']}")
            print(f"Successful: {summary['successful']}")
            print(f"Failed: {summary['failed']}")
            print(f"Success rate: {summary['success_rate']:.1f}%")
            print(f"Output directory: {summary['output_directory']}")
            
            if summary['failed'] > 0:
                return 1
        
        return 0
        
    except Exception as e:
        print(f"Error running simulations: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
