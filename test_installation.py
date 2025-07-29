"""
Test script to verify gensim installation and GCTA availability.
"""

import sys
import os
from pathlib import Path

def test_imports():
    """Test if gensim modules can be imported."""
    try:
        from gensim import GCTASimulator, SimulationConfig, GCTAUtils
        print("✓ Gensim modules imported successfully")
        return True
    except ImportError as e:
        print(f"✗ Failed to import gensim modules: {e}")
        return False

def test_gcta_installation():
    """Test if GCTA is installed and accessible."""
    try:
        from gensim.utils import GCTAUtils
        if GCTAUtils.check_gcta_installation():
            print("✓ GCTA is installed and accessible")
            return True
        else:
            print("✗ GCTA is not accessible (make sure it's in your PATH)")
            return False
    except Exception as e:
        print(f"✗ Error checking GCTA installation: {e}")
        return False

def test_configuration():
    """Test configuration creation."""
    try:
        from gensim import SimulationConfig
        
        config = SimulationConfig(
            bfile="test",
            cohort_sizes=[1000],
            num_causal_snps=[100],
            heritabilities=[0.5]
        )
        
        print("✓ Configuration creation successful")
        print(f"  - Parameter combinations: {len(config.get_parameter_grid())}")
        return True
    except Exception as e:
        print(f"✗ Configuration creation failed: {e}")
        return False

def test_plink_installation():
    """Test PLINK installation."""
    try:
        import subprocess
        result = subprocess.run(
            ["plink", "--version"], 
            capture_output=True, 
            text=True, 
            timeout=10
        )
        if result.returncode == 0:
            print("✓ PLINK found and accessible")
            return True
        else:
            print("⚠ PLINK found but returned error")
            return False
    except Exception as e:
        print(f"✗ PLINK not found or not accessible: {e}")
        print("  Install PLINK from: https://www.cog-genomics.org/plink/")
        return False

def test_plink_simulation():
    """Test PLINK simulation functionality."""
    try:
        from gensim import PLINKSimulator, PLINKSimulationConfig, PLINKSimulationSet
        
        # Create simple test configuration
        snp_sets = [PLINKSimulationSet(100, "test", 0.05, 0.50, 1.00, 1.00)]
        config = PLINKSimulationConfig(
            output_prefix="test_plink",
            num_cases=10,
            num_controls=10,
            snp_sets=snp_sets
        )
        
        print("✓ PLINK simulation configuration successful")
        print(f"  - Total SNPs: {config.get_total_snps()}")
        return True
    except Exception as e:
        print(f"✗ PLINK simulation configuration failed: {e}")
        return False

def main():
    """Run all installation tests."""
    print("Testing gensim installation...")
    print("=" * 40)
    
    results = []
    
    # Test package import
    results.append(test_import())
    
    # Test configuration
    results.append(test_configuration())
    
    # Test GCTA installation
    results.append(test_gcta_installation())
    
    # Test PLINK installation
    results.append(test_plink_installation())
    
    # Test PLINK simulation
    results.append(test_plink_simulation())
    
    # Summary
    print("\n" + "=" * 40)
    passed = sum(results)
    total = len(results)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("✓ All tests passed! Installation looks good.")
        print("\nNext steps:")
        print("1. For GCTA simulations: Prepare your PLINK binary files (.bed, .bim, .fam)")
        print("2. For PLINK dataset creation: Run 'python main.py --create-plink-dataset'")
        print("3. Check examples.py and plink_examples.py for usage examples")
    else:
        print("⚠ Some tests failed. Check the output above for details.")
        if not results[2]:  # GCTA test failed
            print("  - GCTA is required for phenotype simulation")
        if not results[3]:  # PLINK test failed
            print("  - PLINK is required for dataset creation")
    
    return 0 if passed == total else 1
