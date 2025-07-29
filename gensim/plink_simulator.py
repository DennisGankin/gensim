"""
PLINK simulation functionality for generating synthetic SNP datasets.
"""

import os
import subprocess
import logging
from typing import List, Dict, Optional, Tuple, Any, Union
from dataclasses import dataclass, field
from pathlib import Path
from itertools import product

from .utils import GCTAUtils


@dataclass
class PLINKParameterGrid:
    """Configuration for parameter grid generation."""
    
    # Grid parameters - can be single values or lists
    cohort_sizes: Union[int, List[int]] = field(default_factory=lambda: [1000])  # total individuals (cases + controls)
    prevalences: Union[float, List[float]] = field(default_factory=lambda: [0.01])
    total_snps: Union[int, List[int]] = field(default_factory=lambda: [10000])
    causal_snps: Union[int, List[int]] = field(default_factory=lambda: [100])
    
    # Grid output configuration
    grid_output_dir: str = "parameter_grid"
    base_prefix: str = "dataset"
    
    # Fixed simulation parameters (applied to all combinations)
    case_control_ratio: float = 1.0  # 1:1 cases to controls by default
    min_freq: float = 0.001
    max_freq: float = 0.5
    het_odds_ratio: float = 1.5
    hom_odds_ratio: Union[float, str] = "mult"
    plink_executable: str = "plink"
    random_seed: Optional[int] = None
    
    def __post_init__(self):
        """Validate and normalize grid parameters."""
        # Convert single values to lists for consistency
        if isinstance(self.cohort_sizes, int):
            self.cohort_sizes = [self.cohort_sizes]
        if isinstance(self.prevalences, float):
            self.prevalences = [self.prevalences]
        if isinstance(self.total_snps, int):
            self.total_snps = [self.total_snps]
        if isinstance(self.causal_snps, int):
            self.causal_snps = [self.causal_snps]
        
        # Validate parameters
        for cohort_size in self.cohort_sizes:
            if cohort_size < 10:
                raise ValueError(f"Cohort size must be at least 10, got {cohort_size}")
        
        for prevalence in self.prevalences:
            if not 0 < prevalence < 1:
                raise ValueError(f"Prevalence must be between 0 and 1, got {prevalence}")
        
        for total_snp in self.total_snps:
            if total_snp < 1:
                raise ValueError(f"Total SNPs must be positive, got {total_snp}")
        
        for causal_snp in self.causal_snps:
            if causal_snp < 0:
                raise ValueError(f"Causal SNPs cannot be negative, got {causal_snp}")
        
        # Check that causal SNPs don't exceed total SNPs for any combination
        max_causal = max(self.causal_snps)
        min_total = min(self.total_snps)
        if max_causal > min_total:
            raise ValueError(f"Maximum causal SNPs ({max_causal}) exceeds minimum total SNPs ({min_total})")
        
        if not 0 < self.case_control_ratio <= 10:
            raise ValueError(f"Case-control ratio must be between 0 and 10, got {self.case_control_ratio}")
        
        # Create output directory
        os.makedirs(self.grid_output_dir, exist_ok=True)
    
    def get_parameter_combinations(self) -> List[Dict[str, Any]]:
        """
        Generate all parameter combinations from the grid.
        
        Returns:
            List of dictionaries, each containing one parameter combination
        """
        combinations = []
        
        for cohort_size, prevalence, total_snp, causal_snp in product(
            self.cohort_sizes, self.prevalences, self.total_snps, self.causal_snps
        ):
            # Calculate cases and controls based on cohort size and ratio
            total_cases = int(cohort_size / (1 + self.case_control_ratio))
            total_controls = cohort_size - total_cases
            
            # Ensure at least 1 case and 1 control
            if total_cases == 0:
                total_cases = 1
                total_controls = cohort_size - 1
            
            combinations.append({
                'cohort_size': cohort_size,
                'num_cases': total_cases,
                'num_controls': total_controls,
                'prevalence': prevalence,
                'total_snps': total_snp,
                'causal_snps': causal_snp,
                'null_snps': total_snp - causal_snp
            })
        
        return combinations
    
    def get_combination_name(self, combination: Dict[str, Any]) -> str:
        """
        Generate a descriptive name for a parameter combination.
        
        Args:
            combination: Parameter combination dictionary
            
        Returns:
            Descriptive name string
        """
        return (f"{self.base_prefix}_"
                f"n{combination['cohort_size']}_"
                f"prev{combination['prevalence']:.3f}_"
                f"snps{combination['total_snps']}_"
                f"causal{combination['causal_snps']}")


@dataclass
class PLINKSimulationSet:
    """Configuration for a set of SNPs in PLINK simulation."""
    
    num_snps: int
    label: str
    min_freq: float
    max_freq: float
    het_odds_ratio: float
    hom_odds_ratio: float  # or "mult" for multiplicative
    
    def __post_init__(self):
        """Validate simulation set parameters."""
        if not 0 <= self.min_freq <= 1:
            raise ValueError(f"min_freq must be between 0 and 1, got {self.min_freq}")
        if not 0 <= self.max_freq <= 1:
            raise ValueError(f"max_freq must be between 0 and 1, got {self.max_freq}")
        if self.min_freq > self.max_freq:
            raise ValueError(f"min_freq ({self.min_freq}) cannot be greater than max_freq ({self.max_freq})")
        if self.het_odds_ratio < 0:
            raise ValueError(f"het_odds_ratio must be positive, got {self.het_odds_ratio}")
        if isinstance(self.hom_odds_ratio, (int, float)) and self.hom_odds_ratio < 0:
            raise ValueError(f"hom_odds_ratio must be positive, got {self.hom_odds_ratio}")
    
    def to_plink_line(self) -> str:
        """Convert to PLINK simulation file format."""
        hom_or = "mult" if self.hom_odds_ratio == "mult" else f"{self.hom_odds_ratio:.2f}"
        return f"{self.num_snps}\t{self.label}\t{self.min_freq:.2f}\t{self.max_freq:.2f}\t{self.het_odds_ratio:.2f}\t{hom_or}"


@dataclass
class PLINKSimulationConfig:
    """Configuration for PLINK dataset simulation."""
    
    # Output configuration
    output_prefix: str
    output_dir: str = "plink_simulations"
    
    # Population parameters
    num_cases: int = 100
    num_controls: int = 100
    disease_prevalence: float = 0.01
    
    # SNP sets to simulate
    snp_sets: List[PLINKSimulationSet] = field(default_factory=list)
    
    # PLINK executable
    plink_executable: str = "plink"
    
    # Optional parameters
    random_seed: Optional[int] = None
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        if not self.snp_sets:
            raise ValueError("At least one SNP set must be specified")
        
        if not 0 < self.disease_prevalence < 1:
            raise ValueError(f"disease_prevalence must be between 0 and 1, got {self.disease_prevalence}")
        
        if self.num_cases <= 0 or self.num_controls <= 0:
            raise ValueError("num_cases and num_controls must be positive")
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
    
    def get_total_snps(self) -> int:
        """Get total number of SNPs across all sets."""
        return sum(snp_set.num_snps for snp_set in self.snp_sets)
    
    def get_simulation_filename(self) -> str:
        """Get the simulation parameter filename."""
        return f"{self.output_prefix}.sim"
    
    def get_output_path(self, filename: str) -> str:
        """Get full path for output file."""
        return os.path.join(self.output_dir, filename)


class PLINKSimulator:
    """Class for running PLINK simulations to generate synthetic SNP datasets."""
    
    def __init__(self, config: PLINKSimulationConfig):
        """
        Initialize PLINK simulator.
        
        Args:
            config: PLINKSimulationConfig object with simulation parameters
        """
        self.config = config
        self.logger = GCTAUtils.setup_logging(
            log_file=os.path.join(config.output_dir, "plink_simulation.log")
        )
        
        # Validate setup
        self._validate_setup()
    
    def _validate_setup(self):
        """Validate that PLINK is available."""
        if not self._check_plink_installation():
            raise RuntimeError(
                f"PLINK executable '{self.config.plink_executable}' not found. "
                "Please ensure PLINK is installed and accessible."
            )
        
        self.logger.info("PLINK simulation setup validation completed successfully")
    
    def _check_plink_installation(self) -> bool:
        """Check if PLINK is installed and accessible."""
        try:
            result = subprocess.run(
                [self.config.plink_executable, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def create_simulation_file(self) -> str:
        """
        Create PLINK simulation parameter file.
        
        Returns:
            Path to created simulation file
        """
        sim_filename = self.config.get_simulation_filename()
        sim_filepath = self.config.get_output_path(sim_filename)
        
        self.logger.info(f"Creating simulation file: {sim_filepath}")
        
        with open(sim_filepath, 'w') as f:
            for snp_set in self.config.snp_sets:
                f.write(snp_set.to_plink_line() + "\n")
        
        self.logger.info(f"Simulation file created with {len(self.config.snp_sets)} SNP sets")
        return sim_filepath
    
    def run_simulation(self) -> Dict[str, Any]:
        """
        Run PLINK simulation to generate synthetic dataset.
        
        Returns:
            Dictionary with simulation results
        """
        self.logger.info("Starting PLINK simulation")
        
        try:
            # Check if output files already exist and are complete
            output_prefix = self.config.output_prefix
            expected_files = [
                f"{output_prefix}.bed",
                f"{output_prefix}.bim", 
                f"{output_prefix}.fam"
            ]
            
            all_exist = all(
                os.path.exists(os.path.join(self.config.output_dir, f)) 
                for f in expected_files
            )
            
            if all_exist:
                self.logger.info("Output files already exist, skipping simulation")
                
                # Get file statistics for existing simulation
                bed_file = os.path.join(self.config.output_dir, f"{output_prefix}.bed")
                fam_file = os.path.join(self.config.output_dir, f"{output_prefix}.fam")
                bim_file = os.path.join(self.config.output_dir, f"{output_prefix}.bim")
                
                # Count individuals and SNPs from existing files
                try:
                    with open(fam_file, 'r') as f:
                        total_individuals = sum(1 for _ in f)
                    
                    with open(bim_file, 'r') as f:
                        total_snps = sum(1 for _ in f)
                    
                    # Count cases and controls from .fam file (phenotype in column 6)
                    cases = 0
                    controls = 0
                    with open(fam_file, 'r') as f:
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 6:
                                phenotype = parts[5]
                                if phenotype == '2':  # Case
                                    cases += 1
                                elif phenotype == '1':  # Control
                                    controls += 1
                    
                    return {
                        'success': True,
                        'output_prefix': output_prefix,
                        'output_files': expected_files,
                        'skipped': True,
                        'message': 'Simulation already completed',
                        'statistics': {
                            'total_individuals': total_individuals,
                            'total_snps': total_snps,
                            'cases': cases,
                            'controls': controls
                        }
                    }
                
                except Exception as e:
                    self.logger.warning(f"Could not read existing file statistics: {e}")
                    return {
                        'success': True,
                        'output_prefix': output_prefix,
                        'output_files': expected_files,
                        'skipped': True,
                        'message': 'Simulation already completed (statistics unavailable)'
                    }
            
            # Create simulation parameter file
            sim_file = self.create_simulation_file()
            
            # Build PLINK command - use relative paths since we'll run from output_dir
            sim_filename = self.config.get_simulation_filename()
            output_prefix = self.config.output_prefix
            
            cmd = [
                self.config.plink_executable,
                "--simulate", sim_filename,
                "--make-bed",
                "--out", output_prefix,
                "--simulate-ncases", str(self.config.num_cases),
                "--simulate-ncontrols", str(self.config.num_controls),
                "--simulate-prevalence", str(self.config.disease_prevalence)
            ]
            
            # Add random seed if specified
            if self.config.random_seed is not None:
                cmd.extend(["--seed", str(self.config.random_seed)])
            
            # Run PLINK command
            self.logger.info(f"Running command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=14400,  # 4 hour timeout (14400 seconds)
                cwd=self.config.output_dir
            )
            
            success = result.returncode == 0
            
            if success:
                self.logger.info("PLINK simulation completed successfully")
                
                # Check for output files - use full paths
                full_output_prefix = self.config.get_output_path(self.config.output_prefix)
                expected_files = [
                    f"{full_output_prefix}.bed",
                    f"{full_output_prefix}.bim", 
                    f"{full_output_prefix}.fam"
                ]
                
                missing_files = [f for f in expected_files if not os.path.exists(f)]
                
                if missing_files:
                    success = False
                    error_msg = f"Missing expected output files: {missing_files}"
                    self.logger.error(error_msg)
                else:
                    # Get simulation statistics
                    stats = self._get_simulation_stats(full_output_prefix)
                    
                    self.logger.info(f"Generated dataset with:")
                    self.logger.info(f"  - Total SNPs: {stats['total_snps']}")
                    self.logger.info(f"  - Total individuals: {stats['total_individuals']}")
                    self.logger.info(f"  - Cases: {stats['cases']}")
                    self.logger.info(f"  - Controls: {stats['controls']}")
                    
                    return {
                        'success': True,
                        'output_files': expected_files,
                        'output_prefix': full_output_prefix,
                        'simulation_file': sim_file,
                        'statistics': stats,
                        'stdout': result.stdout,
                        'stderr': result.stderr
                    }
            else:
                error_msg = f"PLINK simulation failed with return code {result.returncode}"
                self.logger.error(error_msg)
                self.logger.error(f"STDERR: {result.stderr}")
            
            return {
                'success': success,
                'error': error_msg if not success else None,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'output_files': [],
                'simulation_file': sim_file if 'sim_file' in locals() else None
            }
            
        except Exception as e:
            error_msg = f"Exception during PLINK simulation: {str(e)}"
            self.logger.error(error_msg)
            return {
                'success': False,
                'error': error_msg,
                'output_files': []
            }
    
    def _get_simulation_stats(self, output_prefix: str) -> Dict[str, int]:
        """Get basic statistics about the generated dataset."""
        try:
            import pandas as pd
            
            # Read .fam file for individual counts
            fam_file = f"{output_prefix}.fam"
            fam_df = pd.read_csv(fam_file, sep=r'\s+', header=None,
                               names=['fid', 'iid', 'father', 'mother', 'sex', 'phenotype'])
            
            # Read .bim file for SNP counts
            bim_file = f"{output_prefix}.bim"
            bim_df = pd.read_csv(bim_file, sep='\t', header=None,
                               names=['chr', 'snp_id', 'genetic_dist', 'bp_pos', 'allele1', 'allele2'])
            
            # Count cases and controls (phenotype: 1=control, 2=case)
            cases = len(fam_df[fam_df['phenotype'] == 2])
            controls = len(fam_df[fam_df['phenotype'] == 1])
            
            return {
                'total_snps': len(bim_df),
                'total_individuals': len(fam_df),
                'cases': cases,
                'controls': controls
            }
            
        except Exception as e:
            self.logger.warning(f"Could not get simulation statistics: {e}")
            return {
                'total_snps': -1,
                'total_individuals': -1,
                'cases': -1,
                'controls': -1
            }
    
    @staticmethod
    def create_default_simulation_sets() -> List[PLINKSimulationSet]:
        """Create a set of default SNP simulation parameters."""
        return [
            # Null SNPs with different frequency ranges
            PLINKSimulationSet(20000, "nullA", 0.00, 0.05, 1.00, 1.00),
            PLINKSimulationSet(10000, "nullB", 0.05, 0.10, 1.00, 1.00),
            PLINKSimulationSet(5000, "nullC", 0.10, 0.20, 1.00, 1.00),
            PLINKSimulationSet(10000, "nullD", 0.20, 0.50, 1.00, 1.00),
            # Disease-associated SNPs
            PLINKSimulationSet(100, "disease", 0.05, 0.50, 1.50, "mult"),
        ]
    
    @staticmethod
    def create_simple_simulation_sets(num_null: int = 10000, num_disease: int = 100) -> List[PLINKSimulationSet]:
        """Create simple simulation sets for quick testing."""
        return [
            PLINKSimulationSet(num_null, "null", 0.05, 0.95, 1.00, 1.00),
            PLINKSimulationSet(num_disease, "disease", 0.05, 0.50, 2.00, "mult"),
        ]


class PLINKParameterGridSimulator:
    """Class for running parameter grid simulations with PLINK."""
    
    def __init__(self, grid_config: PLINKParameterGrid):
        """
        Initialize parameter grid simulator.
        
        Args:
            grid_config: PLINKParameterGrid object with grid parameters
        """
        self.grid_config = grid_config
        self.logger = GCTAUtils.setup_logging(
            log_file=os.path.join(grid_config.grid_output_dir, "parameter_grid.log")
        )
        
        # Validate setup
        self._validate_setup()
    
    def _validate_setup(self):
        """Validate that PLINK is available."""
        try:
            result = subprocess.run(
                [self.grid_config.plink_executable, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode != 0:
                raise RuntimeError("PLINK returned non-zero exit code")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            raise RuntimeError(
                f"PLINK executable '{self.grid_config.plink_executable}' not found. "
                "Please ensure PLINK is installed and accessible."
            )
        
        self.logger.info("PLINK parameter grid simulation setup validation completed successfully")
    
    def run_parameter_grid(self) -> Dict[str, Any]:
        """
        Run PLINK simulations for all parameter combinations in the grid.
        
        Returns:
            Dictionary with results for all combinations
        """
        combinations = self.grid_config.get_parameter_combinations()
        total_combinations = len(combinations)
        
        self.logger.info(f"Starting parameter grid simulation with {total_combinations} combinations")
        print(f"Running parameter grid simulation with {total_combinations} combinations...")
        
        results = {
            'total_combinations': total_combinations,
            'successful': 0,
            'failed': 0,
            'combination_results': [],
            'summary': {
                'grid_output_dir': self.grid_config.grid_output_dir,
                'parameters': {
                    'cohort_sizes': self.grid_config.cohort_sizes,
                    'prevalences': self.grid_config.prevalences,
                    'total_snps': self.grid_config.total_snps,
                    'causal_snps': self.grid_config.causal_snps
                }
            }
        }
        
        for i, combination in enumerate(combinations, 1):
            print(f"\nProcessing combination {i}/{total_combinations}:")
            print(f"  Cohort size: {combination['cohort_size']} "
                  f"(Cases: {combination['num_cases']}, Controls: {combination['num_controls']})")
            print(f"  Prevalence: {combination['prevalence']:.3f}")
            print(f"  Total SNPs: {combination['total_snps']} "
                  f"(Causal: {combination['causal_snps']}, Null: {combination['null_snps']})")
            
            try:
                result = self._run_single_combination(combination, i)
                if result['success']:
                    results['successful'] += 1
                    print(f"  ✓ Success: {result['output_prefix']}")
                else:
                    results['failed'] += 1
                    print(f"  ✗ Failed: {result['error']}")
                
                results['combination_results'].append(result)
                
            except Exception as e:
                error_msg = f"Unexpected error in combination {i}: {str(e)}"
                self.logger.error(error_msg)
                print(f"  ✗ Error: {error_msg}")
                
                results['failed'] += 1
                results['combination_results'].append({
                    'combination_number': i,
                    'combination': combination,
                    'success': False,
                    'error': error_msg
                })
        
        # Create summary report
        self._create_summary_report(results)
        
        self.logger.info(f"Parameter grid simulation completed. "
                        f"Success: {results['successful']}, Failed: {results['failed']}")
        
        return results
    
    def _run_single_combination(self, combination: Dict[str, Any], combination_number: int) -> Dict[str, Any]:
        """
        Run PLINK simulation for a single parameter combination.
        
        Args:
            combination: Parameter combination dictionary
            combination_number: Sequential number of this combination
            
        Returns:
            Dictionary with simulation results
        """
        # Generate combination name and output directory
        combination_name = self.grid_config.get_combination_name(combination)
        combination_dir = os.path.join(self.grid_config.grid_output_dir, combination_name)
        
        # Create SNP sets for this combination
        snp_sets = [
            PLINKSimulationSet(
                num_snps=combination['null_snps'],
                label="null",
                min_freq=self.grid_config.min_freq,
                max_freq=self.grid_config.max_freq,
                het_odds_ratio=1.0,
                hom_odds_ratio=1.0
            )
        ]
        
        # Add causal SNPs if any
        if combination['causal_snps'] > 0:
            snp_sets.append(
                PLINKSimulationSet(
                    num_snps=combination['causal_snps'],
                    label="causal",
                    min_freq=self.grid_config.min_freq,
                    max_freq=self.grid_config.max_freq,
                    het_odds_ratio=self.grid_config.het_odds_ratio,
                    hom_odds_ratio=self.grid_config.hom_odds_ratio
                )
            )
        
        # Create PLINK configuration for this combination
        plink_config = PLINKSimulationConfig(
            output_prefix=combination_name,
            output_dir=combination_dir,
            num_cases=combination['num_cases'],
            num_controls=combination['num_controls'],
            disease_prevalence=combination['prevalence'],
            snp_sets=snp_sets,
            plink_executable=self.grid_config.plink_executable,
            random_seed=self.grid_config.random_seed
        )
        
        # Run simulation
        simulator = PLINKSimulator(plink_config)
        simulation_result = simulator.run_simulation()
        
        # Add combination metadata to result
        result = {
            'combination_number': combination_number,
            'combination': combination,
            'combination_name': combination_name,
            'output_directory': combination_dir,
            **simulation_result
        }
        
        return result
    
    def _create_summary_report(self, results: Dict[str, Any]):
        """Create a summary report of the parameter grid simulation."""
        summary_file = os.path.join(self.grid_config.grid_output_dir, "grid_summary.txt")
        
        with open(summary_file, 'w') as f:
            f.write("PLINK Parameter Grid Simulation Summary\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total combinations: {results['total_combinations']}\n")
            f.write(f"Successful: {results['successful']}\n")
            f.write(f"Failed: {results['failed']}\n")
            f.write(f"Success rate: {results['successful']/results['total_combinations']*100:.1f}%\n\n")
            
            f.write("Grid Parameters:\n")
            f.write(f"  Cohort sizes: {results['summary']['parameters']['cohort_sizes']}\n")
            f.write(f"  Prevalences: {results['summary']['parameters']['prevalences']}\n")
            f.write(f"  Total SNPs: {results['summary']['parameters']['total_snps']}\n")
            f.write(f"  Causal SNPs: {results['summary']['parameters']['causal_snps']}\n\n")
            
            f.write("Individual Results:\n")
            f.write("-" * 30 + "\n")
            
            for result in results['combination_results']:
                f.write(f"\nCombination {result['combination_number']}: {result.get('combination_name', 'Unknown')}\n")
                if result['success']:
                    f.write(f"  Status: SUCCESS\n")
                    f.write(f"  Output: {result.get('output_prefix', 'N/A')}\n")
                    if 'statistics' in result:
                        stats = result['statistics']
                        f.write(f"  Statistics: {stats['total_individuals']} individuals, "
                               f"{stats['total_snps']} SNPs, "
                               f"{stats['cases']} cases, {stats['controls']} controls\n")
                else:
                    f.write(f"  Status: FAILED\n")
                    f.write(f"  Error: {result.get('error', 'Unknown error')}\n")
        
        print(f"\nSummary report saved to: {summary_file}")
