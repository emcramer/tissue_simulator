"""
Tissue Simulator Package
========================

A package for creating 3D simulated biological tissue sections using
random sphere packing algorithms, with network-based cell type assignment.
"""

from .tissue import TissueSection, Cell, load_tissue_from_csv
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
    load_target_statistics_from_csv as load_replicate_stats_csv,
    load_target_statistics_from_tissue,
    load_target_statistics_from_coordinates
)
from .graph_coloring import (
    GraphColorizer,
    color_graph_to_targets,
    calculate_graph_statistics,
    compare_graph_statistics,
    load_target_statistics_from_csv as load_graph_stats_csv,
    export_colored_graph_statistics,
    visualize_colored_graph,
    visualize_graph_comparison
)
from .evaluation import (
    cosine_similarity,
    cosine_distance,
    js_divergence,
    evaluate_graph_coloring,
    print_evaluation_report
)
from .tissue_workflow import (
    TissueNetworkWorkflow,
    quick_workflow
)
from .physicell_export import (
    PhysiCellExporter,
    export_to_physicell,
)
from .physicell_reader import (
    PhysiCellReader,
    read_physicell_output,
    stats_to_target_statistics,
)
from .convergence import (
    adf_test,
    mann_kendall_test,
    rolling_cv,
    find_convergence_time,
    MultiMetricConvergence,
)
from .power_analysis import (
    cohens_d,
    coefficient_of_variation,
    required_replicates,
    power_curve,
    compare_initialization_variance,
    summarize_power_analysis,
)

__version__ = "0.1.15"
__all__ = [
    # Core tissue simulation
    "TissueSection", "Cell", "SpherePacker", "load_tissue_from_csv",
    # Slicing
    "TissueSlicer", "SliceCell", "create_standard_slices",
    # Spatial analysis
    "SpatialNetworkAnalyzer", "NetworkStatistics", "CellTypeStatistics",
    "InteractionStatistics", "analyze_tissue_network", "analyze_slice_network",
    # Replicate generation
    "ReplicateGenerator", "TargetStatistics", "ReplicateStatistics",
    "load_replicate_stats_csv", "load_target_statistics_from_tissue",
    "load_target_statistics_from_coordinates",
    # Graph coloring and cell type assignment
    "GraphColorizer", "color_graph_to_targets",
    "calculate_graph_statistics", "compare_graph_statistics",
    "load_graph_stats_csv", "export_colored_graph_statistics",
    "visualize_colored_graph", "visualize_graph_comparison",
    # Evaluation
    "cosine_similarity", "cosine_distance", "js_divergence",
    "evaluate_graph_coloring", "print_evaluation_report",
    # Workflow
    "TissueNetworkWorkflow", "quick_workflow",
    # PhysiCell ABM bridge (v0.1.1)
    "PhysiCellExporter", "export_to_physicell",
    "PhysiCellReader", "read_physicell_output", "stats_to_target_statistics",
    # Convergence diagnostics (v0.1.1)
    "adf_test", "mann_kendall_test", "rolling_cv",
    "find_convergence_time", "MultiMetricConvergence",
    # Power analysis (v0.1.1)
    "cohens_d", "coefficient_of_variation", "required_replicates",
    "power_curve", "compare_initialization_variance",
    "summarize_power_analysis",
]
