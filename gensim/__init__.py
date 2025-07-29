"""
Gensim: A Python package for generating simulated genomic data using GCTA.
"""

__version__ = "0.1.0"
__author__ = "Dennis Gankin"

from .simulator import GCTASimulator
from .config import SimulationConfig
from .utils import GCTAUtils
from .plink_simulator import (
    PLINKSimulator, 
    PLINKSimulationConfig, 
    PLINKSimulationSet,
    PLINKParameterGrid,
    PLINKParameterGridSimulator
)
from .plink_simulator import PLINKSimulator, PLINKSimulationConfig, PLINKSimulationSet

__all__ = ["GCTASimulator", "SimulationConfig", "GCTAUtils", "PLINKSimulator", "PLINKSimulationConfig", "PLINKSimulationSet"]
