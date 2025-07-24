"""
Example script showing how to use gensim for genomic data simulation.
"""

from gensim import GCTASimulator, SimulationConfig


def example_quantitative_trait():
    """Example simulation for quantitative traits."""
    print("Running example quantitative trait simulation...")
    
    # Configure simulation parameters
    config = SimulationConfig(
        bfile="test",  # Base name for your PLINK files
        output_dir="example_quantitative",
        cohort_sizes=[500, 1000, 2000],
        num_causal_snps=[50, 100, 200],
        heritabilities=[0.3, 0.5, 0.8],
        num_replications=3,
        trait_type="quantitative",
        random_seed=42
    )
    
    # Create simulator and run
    simulator = GCTASimulator(config)
    results = simulator.run_simulation_grid()
    
    # Print summary
    summary = simulator.get_results_summary()
    print(f"Completed {summary['successful']}/{summary['total_simulations']} simulations")
    
    return results


def example_binary_trait():
    """Example simulation for binary traits (case-control)."""
    print("Running example binary trait simulation...")
    
    # Configure simulation parameters
    config = SimulationConfig(
        bfile="test",
        output_dir="example_binary",
        cohort_sizes=[1000, 2000],
        num_causal_snps=[100, 200],
        heritabilities=[0.5, 0.8],
        prevalences=[0.05, 0.1, 0.2],  # Disease prevalence
        trait_type="binary",
        random_seed=42
    )
    
    # Create simulator and run
    simulator = GCTASimulator(config)
    results = simulator.run_simulation_grid()
    
    # Print summary
    summary = simulator.get_results_summary()
    print(f"Completed {summary['successful']}/{summary['total_simulations']} simulations")
    
    return results


def example_single_simulation():
    """Example of running a single custom simulation."""
    print("Running single custom simulation...")
    
    # Basic configuration
    config = SimulationConfig(
        bfile="test",
        output_dir="example_single",
        trait_type="quantitative",
        random_seed=42
    )
    
    # Create simulator
    simulator = GCTASimulator(config)
    
    # Run single simulation
    result = simulator.run_single_custom_simulation(
        cohort_size=1000,
        num_causal=100,
        heritability=0.5
    )
    
    if result['success']:
        print(f"Simulation successful! Output files: {result['output_files']}")
    else:
        print(f"Simulation failed: {result['error']}")
    
    return result


def example_with_custom_files():
    """Example using custom causal SNP list and keep files."""
    print("Running simulation with custom files...")
    
    config = SimulationConfig(
        bfile="test",
        output_dir="example_custom",
        cohort_sizes=[1000],
        num_causal_snps=[100],  # This will be ignored if causal_snplist is provided
        heritabilities=[0.5],
        causal_snplist="my_causal_snps.txt",  # Pre-defined causal SNPs
        keep_individuals="my_individuals.keep",  # Pre-defined individuals to keep
        trait_type="quantitative",
        random_seed=42
    )
    
    simulator = GCTASimulator(config)
    results = simulator.run_simulation_grid()
    
    return results


if __name__ == "__main__":
    print("Gensim Example Simulations")
    print("=" * 40)
    
    # Note: Make sure you have PLINK binary files named 'test.bed', 'test.bim', 'test.fam'
    # and GCTA is installed before running these examples
    
    try:
        # Run quantitative trait example
        example_quantitative_trait()
        print()
        
        # Run binary trait example
        example_binary_trait()
        print()
        
        # Run single simulation example
        example_single_simulation()
        print()
        
        print("All examples completed!")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        print("Make sure you have:")
        print("1. GCTA installed and accessible")
        print("2. PLINK binary files (test.bed, test.bim, test.fam)")
        print("3. Proper file permissions")
