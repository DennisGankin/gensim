#!/usr/bin/env python3
"""
Examples of using the parameter grid functionality in Gensim.

This module demonstrates how to create parameter grids for PLINK dataset generation,
allowing systematic exploration of different simulation parameters.
"""

from gensim import PLINKParameterGrid, PLINKParameterGridSimulator


def example_basic_parameter_grid():
    """Basic parameter grid with a few combinations."""
    print("=" * 60)
    print("BASIC PARAMETER GRID EXAMPLE")
    print("=" * 60)
    
    # Define a simple parameter grid
    grid_config = PLINKParameterGrid(
        cohort_sizes=[500, 1000],              # 2 values
        prevalences=[0.01, 0.05],              # 2 values  
        total_snps=[5000, 10000],              # 2 values
        causal_snps=[50, 100],                 # 2 values
        grid_output_dir="basic_grid",
        base_prefix="basic_dataset"
    )
    
    print(f"This will create {2*2*2*2} = 16 dataset combinations")
    
    # Check combinations without running
    combinations = grid_config.get_parameter_combinations()
    print(f"\nFirst few combinations:")
    for i, combo in enumerate(combinations[:3]):
        print(f"  {i+1}. Cohort: {combo['cohort_size']}, "
              f"Prevalence: {combo['prevalence']:.3f}, "
              f"SNPs: {combo['total_snps']}, "
              f"Causal: {combo['causal_snps']}")
    
    print(f"  ... and {len(combinations)-3} more")
    
    # Uncomment to actually run the simulation
    # simulator = PLINKParameterGridSimulator(grid_config)
    # results = simulator.run_parameter_grid()
    # print(f"\nResults: {results['successful']}/{results['total_combinations']} successful")


def example_comprehensive_parameter_grid():
    """Comprehensive parameter grid for thorough analysis."""
    print("=" * 60)
    print("COMPREHENSIVE PARAMETER GRID EXAMPLE")
    print("=" * 60)
    
    # Define a comprehensive parameter grid
    grid_config = PLINKParameterGrid(
        cohort_sizes=[500, 1000, 2000, 5000],        # 4 values
        prevalences=[0.01, 0.05, 0.1, 0.2],         # 4 values
        total_snps=[5000, 10000, 20000, 50000],      # 4 values
        causal_snps=[50, 100, 200, 500],             # 4 values
        grid_output_dir="comprehensive_grid",
        base_prefix="comp_dataset",
        case_control_ratio=1.0,                      # 1:1 case:control ratio
        min_freq=0.01,
        max_freq=0.5,
        het_odds_ratio=1.5,
        hom_odds_ratio="mult",
        random_seed=42
    )
    
    total_combinations = len(grid_config.get_parameter_combinations())
    print(f"This will create {4*4*4*4} = {total_combinations} dataset combinations")
    print("WARNING: This is a large grid and may take significant time to complete!")
    
    # Show some example combinations
    combinations = grid_config.get_parameter_combinations()
    print(f"\nSample combinations:")
    for i in [0, len(combinations)//4, len(combinations)//2, -1]:
        combo = combinations[i]
        name = grid_config.get_combination_name(combo)
        print(f"  {name}")
        print(f"    Cohort: {combo['cohort_size']} individuals "
              f"({combo['num_cases']} cases, {combo['num_controls']} controls)")
        print(f"    Prevalence: {combo['prevalence']:.3f}")
        print(f"    SNPs: {combo['total_snps']} total ({combo['causal_snps']} causal, {combo['null_snps']} null)")
        print()


def example_disease_focused_grid():
    """Parameter grid focused on disease association studies."""
    print("=" * 60)
    print("DISEASE-FOCUSED PARAMETER GRID EXAMPLE")
    print("=" * 60)
    
    grid_config = PLINKParameterGrid(
        cohort_sizes=[1000, 2000, 5000],             # Realistic study sizes
        prevalences=[0.001, 0.01, 0.05, 0.1],       # Range from rare to common diseases
        total_snps=[50000],                          # Fixed large SNP panel
        causal_snps=[10, 50, 100, 200],             # Varying genetic architecture
        grid_output_dir="disease_grid",
        base_prefix="disease_study",
        case_control_ratio=1.0,
        min_freq=0.05,                               # Focus on common variants
        max_freq=0.95,
        het_odds_ratio=2.0,                          # Strong disease effect
        hom_odds_ratio="mult",                       # Multiplicative model
        random_seed=123
    )
    
    combinations = grid_config.get_parameter_combinations()
    print(f"This grid creates {len(combinations)} disease-focused datasets")
    
    # Show power calculation estimates
    print("\nExpected statistical power varies by:")
    print("- Sample size (larger = more power)")
    print("- Disease prevalence (affects case-control balance)")
    print("- Number of causal SNPs (more SNPs = distributed effect)")
    print("- Effect size (odds ratio = 2.0 is strong)")


def example_custom_configuration():
    """Example with custom simulation parameters."""
    print("=" * 60)
    print("CUSTOM CONFIGURATION EXAMPLE")
    print("=" * 60)
    
    grid_config = PLINKParameterGrid(
        cohort_sizes=[800, 1200],
        prevalences=[0.03, 0.07],
        total_snps=[8000, 12000],
        causal_snps=[80, 120],
        grid_output_dir="custom_grid",
        base_prefix="custom_study",
        case_control_ratio=2.0,                      # 2:1 case:control ratio
        min_freq=0.02,                               # Slightly higher MAF threshold
        max_freq=0.4,                                # Lower maximum frequency
        het_odds_ratio=1.3,                          # Moderate effect size
        hom_odds_ratio=1.7,                          # Specific homozygote OR
        plink_executable="plink1.9",                 # Specific PLINK version
        random_seed=456
    )
    
    combinations = grid_config.get_parameter_combinations()
    print(f"Custom grid with {len(combinations)} combinations")
    
    # Show the effect of case_control_ratio
    print(f"\nWith case_control_ratio = 2.0:")
    for combo in combinations[:2]:
        total = combo['cohort_size']
        cases = combo['num_cases']
        controls = combo['num_controls']
        ratio = cases / controls
        print(f"  Cohort {total}: {cases} cases, {controls} controls (ratio: {ratio:.1f})")


def main():
    """Run all parameter grid examples."""
    print("GENSIM Parameter Grid Examples")
    print("=" * 80)
    print("This script demonstrates different ways to set up parameter grids")
    print("for systematic PLINK dataset generation.\n")
    
    example_basic_parameter_grid()
    print("\n")
    
    example_comprehensive_parameter_grid()
    print("\n")
    
    example_disease_focused_grid()
    print("\n")
    
    example_custom_configuration()
    
    print("\n" + "=" * 80)
    print("To run any of these examples:")
    print("1. Uncomment the simulator.run_parameter_grid() lines")
    print("2. Ensure PLINK is installed and accessible")
    print("3. Run this script or use the CLI interface")
    print("\nCLI example:")
    print("python main.py --create-parameter-grid \\")
    print("  --grid-cohort-sizes '500,1000' \\")
    print("  --grid-prevalences '0.01,0.05' \\")
    print("  --grid-total-snps '5000,10000' \\")
    print("  --grid-causal-snps '50,100'")


if __name__ == "__main__":
    main()
