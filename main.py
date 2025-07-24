#!/usr/bin/env python3
"""
Command-line interface for gensim - genomic data simulation using GCTA.
"""

import argparse
import json
import sys
from pathlib import Path

from gensim import GCTASimulator, SimulationConfig


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
    
    args = parser.parse_args()
    
    # Handle create-config option
    if args.create_config:
        create_example_config(args.create_config)
        return 0
    
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
            print("Error: --bfile is required when not using --config")
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
