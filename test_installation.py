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

def test_sample_files():
    """Check if sample PLINK files exist."""
    bfile = "test"
    extensions = [".bed", ".bim", ".fam"]
    
    missing_files = []
    for ext in extensions:
        if not os.path.exists(f"{bfile}{ext}"):
            missing_files.append(f"{bfile}{ext}")
    
    if missing_files:
        print(f"⚠ Sample PLINK files not found: {missing_files}")
        print("  This is expected if you haven't set up test data yet")
        return False
    else:
        print("✓ Sample PLINK files found")
        return True

def main():
    """Run all tests."""
    print("Gensim Installation Test")
    print("=" * 30)
    
    tests = [
        ("Import Test", test_imports),
        ("GCTA Installation", test_gcta_installation),
        ("Configuration Test", test_configuration),
        ("Sample Files", test_sample_files)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        if test_func():
            passed += 1
    
    print(f"\n{'=' * 30}")
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("✓ All tests passed! Gensim is ready to use.")
        return 0
    else:
        print("⚠ Some tests failed. Check the output above for details.")
        
        if passed >= 2:  # Import and config tests passed
            print("\nYou can still use gensim, but you may need to:")
            print("- Install GCTA and add it to your PATH")
            print("- Provide PLINK binary files for simulation")
        
        return 1

if __name__ == "__main__":
    sys.exit(main())
