"""
Utility functions for GCTA simulation workflows.
"""

import os
import random
import subprocess
import logging
import pandas as pd
from typing import List, Optional, Tuple
from pathlib import Path


class GCTAUtils:
    """Utility class for GCTA-related operations."""
    
    @staticmethod
    def setup_logging(log_file: Optional[str] = None, level: str = "INFO") -> logging.Logger:
        """Set up logging for simulation runs."""
        logger = logging.getLogger("gensim")
        logger.setLevel(getattr(logging, level.upper()))
        
        # Clear existing handlers
        logger.handlers.clear()
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        # File handler if specified
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        
        return logger
    
    @staticmethod
    def check_gcta_installation(executable: str = "gcta64") -> bool:
        """Check if GCTA is installed and accessible."""
        try:
            result = subprocess.run(
                [executable, "--help"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    @staticmethod
    def check_bfile_exists(bfile: str) -> Tuple[bool, List[str]]:
        """
        Check if PLINK binary files exist (.bed, .bim, .fam).
        
        Returns:
            Tuple of (all_exist, missing_files)
        """
        extensions = [".bed", ".bim", ".fam"]
        missing_files = []
        
        for ext in extensions:
            filepath = f"{bfile}{ext}"
            if not os.path.exists(filepath):
                missing_files.append(filepath)
        
        return len(missing_files) == 0, missing_files
    
    @staticmethod
    def create_causal_snplist(bim_file: str, num_causal: int, 
                            output_file: str, random_seed: Optional[int] = None) -> str:
        """
        Create a file with randomly selected causal SNPs from a .bim file.
        
        Args:
            bim_file: Path to .bim file
            num_causal: Number of causal SNPs to select
            output_file: Output file path for SNP list
            random_seed: Random seed for reproducibility
            
        Returns:
            Path to created SNP list file
        """
        if random_seed is not None:
            random.seed(random_seed)
        
        # Read SNP IDs from .bim file
        try:
            bim_df = pd.read_csv(
                bim_file, 
                sep='\t', 
                header=None,
                names=['chr', 'snp_id', 'genetic_dist', 'bp_pos', 'allele1', 'allele2']
            )
        except Exception as e:
            raise ValueError(f"Error reading .bim file {bim_file}: {e}")
        
        if len(bim_df) < num_causal:
            raise ValueError(
                f"Requested {num_causal} causal SNPs but only {len(bim_df)} SNPs available"
            )
        
        # Randomly select causal SNPs
        causal_snps = random.sample(list(bim_df['snp_id']), num_causal)
        
        # Create output directory if needed
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write SNP list
        with open(output_file, 'w') as f:
            for snp in causal_snps:
                f.write(f"{snp}\n")
        
        return output_file
    
    @staticmethod
    def create_keep_file(fam_file: str, cohort_size: int, 
                        output_file: str, random_seed: Optional[int] = None) -> str:
        """
        Create a file with individuals to keep for simulation.
        
        Args:
            fam_file: Path to .fam file
            cohort_size: Number of individuals to keep
            output_file: Output file path for keep list
            random_seed: Random seed for reproducibility
            
        Returns:
            Path to created keep file
        """
        if random_seed is not None:
            random.seed(random_seed)
        
        # Read individuals from .fam file
        try:
            fam_df = pd.read_csv(
                fam_file, 
                sep='\s+', 
                header=None,
                names=['fid', 'iid', 'father', 'mother', 'sex', 'phenotype']
            )
        except Exception as e:
            raise ValueError(f"Error reading .fam file {fam_file}: {e}")
        
        if len(fam_df) < cohort_size:
            raise ValueError(
                f"Requested cohort size {cohort_size} but only {len(fam_df)} individuals available"
            )
        
        # Randomly select individuals
        selected_indices = random.sample(range(len(fam_df)), cohort_size)
        selected_individuals = fam_df.iloc[selected_indices]
        
        # Create output directory if needed
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write keep file (FID IID format)
        selected_individuals[['fid', 'iid']].to_csv(
            output_file, 
            sep='\t', 
            header=False, 
            index=False
        )
        
        return output_file
    
    @staticmethod
    def run_command(command: List[str], logger: Optional[logging.Logger] = None) -> Tuple[bool, str, str]:
        """
        Run a command and return success status and output.
        
        Args:
            command: Command as list of strings
            logger: Optional logger for output
            
        Returns:
            Tuple of (success, stdout, stderr)
        """
        if logger:
            logger.info(f"Running command: {' '.join(command)}")
        
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            success = result.returncode == 0
            
            if logger:
                if success:
                    logger.info("Command completed successfully")
                else:
                    logger.error(f"Command failed with return code {result.returncode}")
                
                if result.stdout:
                    logger.info(f"STDOUT: {result.stdout}")
                if result.stderr:
                    logger.warning(f"STDERR: {result.stderr}")
            
            return success, result.stdout, result.stderr
            
        except subprocess.TimeoutExpired:
            error_msg = "Command timed out after 1 hour"
            if logger:
                logger.error(error_msg)
            return False, "", error_msg
        except Exception as e:
            error_msg = f"Error running command: {e}"
            if logger:
                logger.error(error_msg)
            return False, "", error_msg
    
    @staticmethod
    def validate_simulation_outputs(output_prefix: str, trait_type: str) -> List[str]:
        """
        Validate that expected simulation output files were created.
        
        Args:
            output_prefix: Output prefix used in GCTA command
            trait_type: Type of trait simulated
            
        Returns:
            List of missing files (empty if all files exist)
        """
        expected_files = []
        
        if trait_type == "quantitative":
            expected_files.extend([
                f"{output_prefix}.par",
                f"{output_prefix}.phen"
            ])
        elif trait_type == "binary":
            expected_files.extend([
                f"{output_prefix}.par",
                f"{output_prefix}.phen"
            ])
        
        # Check for log file
        expected_files.append(f"{output_prefix}.log")
        
        missing_files = []
        for file_path in expected_files:
            if not os.path.exists(file_path):
                missing_files.append(file_path)
        
        return missing_files
    
    @staticmethod
    def cleanup_temp_files(file_paths: List[str], logger: Optional[logging.Logger] = None):
        """Remove temporary files."""
        for file_path in file_paths:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    if logger:
                        logger.debug(f"Removed temporary file: {file_path}")
            except Exception as e:
                if logger:
                    logger.warning(f"Failed to remove {file_path}: {e}")
    
    @staticmethod
    def create_summary_report(output_dir: str, simulations: List[dict]) -> str:
        """
        Create a summary report of all simulations.
        
        Args:
            output_dir: Directory containing simulation outputs
            simulations: List of simulation parameters and results
            
        Returns:
            Path to summary report file
        """
        report_file = os.path.join(output_dir, "simulation_summary.csv")
        
        # Convert simulations to DataFrame
        df = pd.DataFrame(simulations)
        
        # Save to CSV
        df.to_csv(report_file, index=False)
        
        return report_file
