"""
Main GCTA simulation class for generating genomic data with parameter grids.
"""

import os
import time
from typing import List, Dict, Optional, Any
from datetime import datetime

from .config import SimulationConfig, GCTACommand
from .utils import GCTAUtils


class GCTASimulator:
    """Main class for running GCTA simulations with parameter grids."""
    
    def __init__(self, config: SimulationConfig):
        """
        Initialize the GCTA simulator.
        
        Args:
            config: SimulationConfig object with all parameters
        """
        self.config = config
        self.logger = GCTAUtils.setup_logging(
            log_file=os.path.join(config.output_dir, "gensim.log")
        )
        self.simulation_results = []
        
        # Validate setup
        self._validate_setup()
    
    def _validate_setup(self):
        """Validate that all required components are available."""
        # Check GCTA installation
        if not GCTAUtils.check_gcta_installation(self.config.gcta_executable):
            raise RuntimeError(
                f"GCTA executable '{self.config.gcta_executable}' not found. "
                "Please ensure GCTA is installed and accessible."
            )
        
        # Check PLINK binary files
        files_exist, missing_files = GCTAUtils.check_bfile_exists(self.config.bfile)
        if not files_exist:
            raise FileNotFoundError(
                f"Missing PLINK binary files: {missing_files}"
            )
        
        self.logger.info("Setup validation completed successfully")
    
    def run_simulation_grid(self, cleanup_temp: bool = True) -> List[Dict[str, Any]]:
        """
        Run simulations for all parameter combinations in the grid.
        
        Args:
            cleanup_temp: Whether to clean up temporary files after simulation
            
        Returns:
            List of dictionaries containing simulation results
        """
        self.logger.info("Starting simulation grid")
        parameter_grid = self.config.get_parameter_grid()
        
        self.logger.info(f"Total parameter combinations: {len(parameter_grid)}")
        
        start_time = time.time()
        successful_sims = 0
        failed_sims = 0
        
        for i, params in enumerate(parameter_grid, 1):
            self.logger.info(
                f"Running simulation {i}/{len(parameter_grid)} - "
                f"Cohort: {params['cohort_size']}, "
                f"Causal SNPs: {params['num_causal']}, "
                f"Heritability: {params['heritability']}"
                + (f", Prevalence: {params['prevalence']}" if params['prevalence'] else "")
            )
            
            try:
                result = self._run_single_simulation(params, cleanup_temp)
                self.simulation_results.append(result)
                
                if result['success']:
                    successful_sims += 1
                    self.logger.info(f"Simulation {i} completed successfully")
                else:
                    failed_sims += 1
                    self.logger.error(f"Simulation {i} failed: {result['error']}")
                    
            except Exception as e:
                failed_sims += 1
                error_result = {
                    **params,
                    'success': False,
                    'error': str(e),
                    'output_files': [],
                    'simulation_name': 'failed',
                    'timestamp': datetime.now().isoformat()
                }
                self.simulation_results.append(error_result)
                self.logger.error(f"Simulation {i} failed with exception: {e}")
        
        total_time = time.time() - start_time
        
        self.logger.info(
            f"Simulation grid completed in {total_time:.2f} seconds. "
            f"Successful: {successful_sims}, Failed: {failed_sims}"
        )
        
        # Create summary report
        summary_file = GCTAUtils.create_summary_report(
            self.config.output_dir, 
            self.simulation_results
        )
        self.logger.info(f"Summary report saved to: {summary_file}")
        
        return self.simulation_results
    
    def _run_single_simulation(self, params: Dict[str, Any], 
                             cleanup_temp: bool = True) -> Dict[str, Any]:
        """
        Run a single simulation with given parameters.
        
        Args:
            params: Dictionary with simulation parameters
            cleanup_temp: Whether to clean up temporary files
            
        Returns:
            Dictionary with simulation results
        """
        temp_files = []
        
        try:
            # Generate simulation name
            sim_name = self.config.get_simulation_name(
                cohort_size=params['cohort_size'],
                num_causal=params['num_causal'],
                heritability=params['heritability'],
                prevalence=params['prevalence']
            )
            
            output_prefix = os.path.join(self.config.output_dir, sim_name)
            
            # Create causal SNP list if not provided
            causal_snplist = self.config.causal_snplist
            if causal_snplist is None:
                causal_snplist = f"{output_prefix}.causal.snplist"
                GCTAUtils.create_causal_snplist(
                    bim_file=f"{self.config.bfile}.bim",
                    num_causal=params['num_causal'],
                    output_file=causal_snplist,
                    random_seed=self.config.random_seed
                )
                temp_files.append(causal_snplist)
            
            # Create keep file for cohort size if needed
            keep_file = self.config.keep_individuals
            if params['cohort_size'] != self._get_total_individuals():
                keep_file = f"{output_prefix}.keep"
                GCTAUtils.create_keep_file(
                    fam_file=f"{self.config.bfile}.fam",
                    cohort_size=params['cohort_size'],
                    output_file=keep_file,
                    random_seed=self.config.random_seed
                )
                temp_files.append(keep_file)
            
            # Build GCTA command
            gcta_cmd = GCTACommand(
                executable=self.config.gcta_executable,
                bfile=self.config.bfile,
                output=output_prefix,
                trait_type=self.config.trait_type,
                heritability=params['heritability'],
                num_replications=self.config.num_replications,
                causal_snplist=causal_snplist,
                prevalence=params['prevalence'],
                keep_individuals=keep_file,
                random_seed=self.config.random_seed
            )
            
            # Run GCTA command
            success, stdout, stderr = GCTAUtils.run_command(
                gcta_cmd.build_command(),
                logger=self.logger
            )
            
            # Validate outputs
            missing_files = GCTAUtils.validate_simulation_outputs(
                output_prefix, 
                self.config.trait_type
            )
            
            if missing_files:
                success = False
                error_msg = f"Missing output files: {missing_files}"
            else:
                error_msg = None
            
            # Collect output files
            output_files = []
            for ext in ['.par', '.phen', '.log']:
                file_path = f"{output_prefix}{ext}"
                if os.path.exists(file_path):
                    output_files.append(file_path)
            
            result = {
                **params,
                'success': success,
                'error': error_msg or stderr if not success else None,
                'output_files': output_files,
                'simulation_name': sim_name,
                'gcta_command': gcta_cmd.command_string(),
                'timestamp': datetime.now().isoformat()
            }
            
            return result
            
        finally:
            # Cleanup temporary files if requested
            if cleanup_temp and temp_files:
                GCTAUtils.cleanup_temp_files(temp_files, self.logger)
    
    def _get_total_individuals(self) -> int:
        """Get total number of individuals in the dataset."""
        try:
            import pandas as pd
            fam_df = pd.read_csv(
                f"{self.config.bfile}.fam", 
                sep='\s+', 
                header=None
            )
            return len(fam_df)
        except Exception:
            # If we can't read the file, assume we need all individuals
            return float('inf')
    
    def get_results_summary(self) -> Dict[str, Any]:
        """Get a summary of simulation results."""
        if not self.simulation_results:
            return {"message": "No simulations have been run yet"}
        
        successful = [r for r in self.simulation_results if r['success']]
        failed = [r for r in self.simulation_results if not r['success']]
        
        return {
            "total_simulations": len(self.simulation_results),
            "successful": len(successful),
            "failed": len(failed),
            "success_rate": len(successful) / len(self.simulation_results) * 100,
            "parameter_combinations": len(self.config.get_parameter_grid()),
            "output_directory": self.config.output_dir
        }
    
    def run_single_custom_simulation(self, cohort_size: int, num_causal: int,
                                   heritability: float, prevalence: Optional[float] = None,
                                   output_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Run a single simulation with custom parameters.
        
        Args:
            cohort_size: Number of individuals
            num_causal: Number of causal SNPs
            heritability: Heritability value
            prevalence: Prevalence for binary traits
            output_name: Custom output name (optional)
            
        Returns:
            Dictionary with simulation results
        """
        params = {
            'cohort_size': cohort_size,
            'num_causal': num_causal,
            'heritability': heritability,
            'prevalence': prevalence
        }
        
        self.logger.info(f"Running custom simulation with parameters: {params}")
        
        result = self._run_single_simulation(params)
        
        # Optionally rename output files
        if output_name and result['success']:
            # Implementation for renaming files would go here
            pass
        
        return result
