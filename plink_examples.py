"""
Example script showing how to use PLINK simulation functionality in gensim.
"""

from gensim import PLINKSimulator, PLINKSimulationConfig, PLINKSimulationSet


def example_simple_plink_dataset():
    """Example of creating a simple PLINK dataset."""
    print("Creating simple PLINK dataset...")
    
    # Create SNP sets
    snp_sets = PLINKSimulator.create_simple_simulation_sets(
        num_null=5000,    # 5,000 null SNPs
        num_disease=50    # 50 disease-associated SNPs
    )
    
    # Configure simulation
    config = PLINKSimulationConfig(
        output_prefix="simple_dataset",
        output_dir="plink_examples",
        num_cases=500,
        num_controls=500,
        disease_prevalence=0.05,
        snp_sets=snp_sets,
        random_seed=42
    )
    
    # Run simulation
    simulator = PLINKSimulator(config)
    result = simulator.run_simulation()
    
    if result['success']:
        print(f"✓ Dataset created: {result['output_prefix']}")
        print(f"  Files: {result['output_files']}")
        stats = result.get('statistics', {})
        print(f"  SNPs: {stats.get('total_snps', 'N/A')}")
        print(f"  Individuals: {stats.get('total_individuals', 'N/A')}")
    else:
        print(f"✗ Failed: {result['error']}")
    
    return result


def example_custom_plink_dataset():
    """Example of creating a custom PLINK dataset with multiple SNP sets."""
    print("Creating custom PLINK dataset with multiple SNP frequency ranges...")
    
    # Create custom SNP sets with different frequency ranges
    snp_sets = [
        # Rare variants (MAF < 5%)
        PLINKSimulationSet(10000, "rare_null", 0.001, 0.05, 1.00, 1.00),
        PLINKSimulationSet(20, "rare_disease", 0.001, 0.05, 2.50, "mult"),
        
        # Common variants (MAF 5-50%)
        PLINKSimulationSet(15000, "common_null", 0.05, 0.50, 1.00, 1.00),
        PLINKSimulationSet(80, "common_disease", 0.05, 0.50, 1.80, "mult"),
        
        # Very common variants (MAF > 20%)
        PLINKSimulationSet(5000, "frequent_null", 0.20, 0.50, 1.00, 1.00),
        PLINKSimulationSet(30, "frequent_disease", 0.20, 0.50, 1.50, 2.25),  # Specific homozygote OR
    ]
    
    # Configure simulation
    config = PLINKSimulationConfig(
        output_prefix="custom_dataset",
        output_dir="plink_examples",
        num_cases=1000,
        num_controls=2000,
        disease_prevalence=0.01,
        snp_sets=snp_sets,
        random_seed=123
    )
    
    # Run simulation
    simulator = PLINKSimulator(config)
    result = simulator.run_simulation()
    
    if result['success']:
        print(f"✓ Custom dataset created: {result['output_prefix']}")
        print(f"  Total SNP sets: {len(snp_sets)}")
        print(f"  Total SNPs: {config.get_total_snps()}")
        stats = result.get('statistics', {})
        print(f"  Actual SNPs in output: {stats.get('total_snps', 'N/A')}")
        print(f"  Cases: {stats.get('cases', 'N/A')}")
        print(f"  Controls: {stats.get('controls', 'N/A')}")
    else:
        print(f"✗ Failed: {result['error']}")
    
    return result


def example_realistic_gwas_dataset():
    """Example of creating a more realistic GWAS dataset."""
    print("Creating realistic GWAS dataset...")
    
    # Use default simulation sets (more realistic frequency distribution)
    snp_sets = PLINKSimulator.create_default_simulation_sets()
    
    # Configure for typical GWAS
    config = PLINKSimulationConfig(
        output_prefix="gwas_dataset",
        output_dir="plink_examples",
        num_cases=2000,
        num_controls=3000,
        disease_prevalence=0.05,  # 5% disease prevalence
        snp_sets=snp_sets,
        random_seed=456
    )
    
    # Run simulation
    simulator = PLINKSimulator(config)
    result = simulator.run_simulation()
    
    if result['success']:
        print(f"✓ GWAS dataset created: {result['output_prefix']}")
        print(f"  Total SNPs: {config.get_total_snps()}")
        print(f"  SNP sets: {len(snp_sets)}")
        for snp_set in snp_sets:
            print(f"    - {snp_set.label}: {snp_set.num_snps} SNPs (MAF: {snp_set.min_freq}-{snp_set.max_freq})")
    else:
        print(f"✗ Failed: {result['error']}")
    
    return result


if __name__ == "__main__":
    print("PLINK Simulation Examples")
    print("=" * 40)
    
    # Note: Make sure PLINK is installed before running these examples
    
    try:
        # Run simple example
        print("\n1. Simple Dataset:")
        example_simple_plink_dataset()
        
        print("\n2. Custom Dataset:")
        example_custom_plink_dataset()
        
        print("\n3. Realistic GWAS Dataset:")
        example_realistic_gwas_dataset()
        
        print("\n✓ All examples completed!")
        
    except Exception as e:
        print(f"✗ Error running examples: {e}")
        print("Make sure PLINK is installed and accessible.")
