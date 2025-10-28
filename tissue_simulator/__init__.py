"""
Tissue Simulator Package
========================

A package for creating 3D simulated biological tissue sections using
random sphere packing algorithms.
"""

from .tissue import TissueSection, Cell
from .packing import SpherePacker
from .slicing import TissueSlicer, SliceCell, create_standard_slices
from .spatial_analysis import (
    SpatialNetworkAnalyzer,
    NetworkStatistics,
    CellTypeStatistics,
    InteractionStatistics,
    analyze_tissue_network,
    analyze_slice_network
)

__version__ = "0.1.0"
__all__ = [
    "TissueSection", "Cell", "SpherePacker",
    "TissueSlicer", "SliceCell", "create_standard_slices",
    "SpatialNetworkAnalyzer", "NetworkStatistics", "CellTypeStatistics",
    "InteractionStatistics", "analyze_tissue_network", "analyze_slice_network"
]
