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

# Import HDF5 utilities if h5py is available
try:
    from .h5_utils import H5PLINKReader, read_h5_plink, list_h5_files
    HDF5_AVAILABLE = True
    __all__ = [
        "GCTASimulator", "SimulationConfig", "GCTAUtils", 
        "PLINKSimulator", "PLINKSimulationConfig", "PLINKSimulationSet",
        "PLINKParameterGrid", "PLINKParameterGridSimulator",
        "H5PLINKReader", "read_h5_plink", "list_h5_files"
    ]
except ImportError:
    HDF5_AVAILABLE = False
    __all__ = [
        "GCTASimulator", "SimulationConfig", "GCTAUtils", 
        "PLINKSimulator", "PLINKSimulationConfig", "PLINKSimulationSet",
        "PLINKParameterGrid", "PLINKParameterGridSimulator"
    ]
