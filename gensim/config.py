"""
Configuration classes for GCTA simulation parameters.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Union
import itertools
import os


@dataclass
class SimulationConfig:
    """Configuration class for GCTA simulation parameters."""
    
    # Basic required parameters
    bfile: str  # Base name for PLINK binary files (.bed, .bim, .fam)
    output_dir: str = "simulations"
    
    # Simulation parameters with grid support
    cohort_sizes: List[int] = field(default_factory=lambda: [1000])
    num_causal_snps: List[int] = field(default_factory=lambda: [100])
    heritabilities: List[float] = field(default_factory=lambda: [0.5])
    prevalences: List[float] = field(default_factory=lambda: [0.1])  # For binary traits
    
    # GCTA-specific parameters
    num_replications: int = 1
    trait_type: str = "quantitative"  # "quantitative" or "binary"
    gcta_executable: str = "gcta64"
    
    # Optional files
    causal_snplist: Optional[str] = None  # If None, will be generated
    keep_individuals: Optional[str] = None  # File with individuals to keep
    
    # Random seed for reproducibility
    random_seed: Optional[int] = None
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        if self.trait_type not in ["quantitative", "binary"]:
            raise ValueError("trait_type must be 'quantitative' or 'binary'")
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Validate heritability values
        for h in self.heritabilities:
            if not 0 <= h <= 1:
                raise ValueError(f"Heritability must be between 0 and 1, got {h}")
        
        # Validate prevalence values for binary traits
        if self.trait_type == "binary":
            for p in self.prevalences:
                if not 0 < p < 1:
                    raise ValueError(f"Prevalence must be between 0 and 1, got {p}")
    
    def get_parameter_grid(self):
        """Generate all combinations of simulation parameters."""
        if self.trait_type == "quantitative":
            # For quantitative traits, we don't use prevalence
            combinations = list(itertools.product(
                self.cohort_sizes,
                self.num_causal_snps,
                self.heritabilities,
                [None]  # Placeholder for prevalence
            ))
        else:
            # For binary traits, include prevalence
            combinations = list(itertools.product(
                self.cohort_sizes,
                self.num_causal_snps,
                self.heritabilities,
                self.prevalences
            ))
        
        return [
            {
                "cohort_size": combo[0],
                "num_causal": combo[1],
                "heritability": combo[2],
                "prevalence": combo[3]
            }
            for combo in combinations
        ]
    
    def get_simulation_name(self, cohort_size: int, num_causal: int, 
                          heritability: float, prevalence: Optional[float] = None,
                          rep: int = 1) -> str:
        """Generate a descriptive name for a simulation run."""
        base_name = f"sim_n{cohort_size}_causal{num_causal}_h{heritability:.2f}"
        
        if self.trait_type == "binary" and prevalence is not None:
            base_name += f"_prev{prevalence:.3f}"
        
        if self.num_replications > 1:
            base_name += f"_rep{rep}"
        
        return base_name


@dataclass
class GCTACommand:
    """Class to build and store GCTA command parameters."""
    
    executable: str
    bfile: str
    output: str
    trait_type: str
    heritability: float
    num_replications: int = 1
    causal_snplist: Optional[str] = None
    prevalence: Optional[float] = None
    keep_individuals: Optional[str] = None
    random_seed: Optional[int] = None
    
    def build_command(self) -> List[str]:
        """Build the complete GCTA command as a list of arguments."""
        cmd = [
            self.executable,
            "--bfile", self.bfile,
            "--out", self.output
        ]
        
        # Add trait-specific parameters
        if self.trait_type == "quantitative":
            cmd.extend(["--simu-qt"])
        elif self.trait_type == "binary":
            cmd.extend(["--simu-cc"])
            if self.prevalence is not None:
                cmd.extend(["--simu-prevalence", str(self.prevalence)])
        
        # Add heritability
        cmd.extend(["--simu-hsq", str(self.heritability)])
        
        # Add causal SNPs
        if self.causal_snplist:
            cmd.extend(["--simu-causal-loci", self.causal_snplist])
        
        # Add number of replications
        if self.num_replications > 1:
            cmd.extend(["--simu-rep", str(self.num_replications)])
        
        # Add individuals to keep
        if self.keep_individuals:
            cmd.extend(["--keep", self.keep_individuals])
        
        # Add random seed
        if self.random_seed is not None:
            cmd.extend(["--seed", str(self.random_seed)])
        
        return cmd
    
    def command_string(self) -> str:
        """Return the command as a single string."""
        return " ".join(self.build_command())
