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
from .replicate_generator import (
    ReplicateGenerator,
    TargetStatistics,
    ReplicateStatistics,
    load_target_statistics_from_csv,
    load_target_statistics_from_tissue
)

__version__ = "0.1.0"
__all__ = [
    "TissueSection", "Cell", "SpherePacker",
    "TissueSlicer", "SliceCell", "create_standard_slices",
    "SpatialNetworkAnalyzer", "NetworkStatistics", "CellTypeStatistics",
    "InteractionStatistics", "analyze_tissue_network", "analyze_slice_network",
    "ReplicateGenerator", "TargetStatistics", "ReplicateStatistics",
    "load_target_statistics_from_csv", "load_target_statistics_from_tissue"
]
