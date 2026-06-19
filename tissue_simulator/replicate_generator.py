"""
Replicate tissue generator based on spatial statistics.

This module allows generation of random tissue replicates that match specified
spatial interaction statistics. It uses an optimization-based approach to tune
tissue parameters to achieve desired cell-cell interaction patterns.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
import csv
from collections import defaultdict
from dataclasses import dataclass, asdict
import warnings
from pathlib import Path

from .tissue import TissueSection, Cell, load_tissue_from_csv
from .packing import SpherePacker
from .spatial_analysis import SpatialNetworkAnalyzer, InteractionStatistics
from .graph_coloring import GraphColorizer, color_graph_to_targets
from .power_analysis import compare_initialization_variance

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    warnings.warn("NetworkX not installed. Install with: pip install networkx")


@dataclass
class TargetStatistics:
    """
    Target spatial statistics for replicate generation.
    
    Attributes:
        interaction_stats: List of InteractionStatistics defining target patterns
        cell_type_proportions: Dict mapping cell types to target proportions (0-1)
        target_cell_count: Optional target for total cell count
        target_density: Optional target for packing fraction
    """
    interaction_stats: List[InteractionStatistics]
    cell_type_proportions: Optional[Dict[str, float]] = None
    target_cell_count: Optional[int] = None
    target_density: Optional[float] = None
    
    def validate(self):
        """Validate that statistics are consistent."""
        if self.cell_type_proportions:
            total_prop = sum(self.cell_type_proportions.values())
            if not (0.99 <= total_prop <= 1.01):
                raise ValueError(f"Cell type proportions must sum to 1.0, got {total_prop}")
        
        if self.target_density and not (0 < self.target_density < 1):
            raise ValueError(f"Target density must be between 0 and 1, got {self.target_density}")


@dataclass
class ReplicateStatistics:
    """
    Statistics for a generated replicate.
    
    Attributes:
        replicate_id: Unique identifier
        num_cells: Total cell count
        cell_type_counts: Dict of counts per cell type
        packing_fraction: Volume fraction occupied by cells
        interaction_stats: Measured interaction statistics
        divergence_score: Overall divergence from target statistics
    """
    replicate_id: int
    num_cells: int
    cell_type_counts: Dict[str, int]
    packing_fraction: float
    interaction_stats: List[InteractionStatistics]
    divergence_score: float
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        result = asdict(self)
        result['interaction_stats'] = [s.to_dict() for s in self.interaction_stats]
        return result


class ReplicateGenerator:
    """
    Generate tissue replicates matching target spatial statistics.
    
    This class uses an iterative approach to generate tissue samples that
    match specified spatial interaction patterns. It adjusts tissue parameters
    to achieve the desired statistics.
    """
    
    def __init__(self, 
                 target_stats: TargetStatistics,
                 tissue_dimensions: Tuple[float, float, float],
                 base_cell_radii: Dict[str, Tuple[float, float]],
                 network_mode: str = "contact",
                 network_radius: Optional[float] = None,
                 seed: Optional[int] = None,
                 method: str = "radius_tuning",
                 coloring_params: Optional[Dict] = None,
                 n_restarts: int = 1,
                 radius_optimizer: str = "heuristic",
                 de_params: Optional[Dict] = None):
        """
        Initialize replicate generator.

        Args:
            target_stats: Target spatial statistics to match
            tissue_dimensions: (height, width, thickness) in micrometers
            base_cell_radii: Dict mapping cell types to (min_radius, max_radius)
            network_mode: "contact" or "radius" for spatial analysis
            network_radius: Distance threshold if using "radius" mode
            seed: Random seed for reproducibility
            method: Replicate strategy. ``"radius_tuning"`` (default, unchanged
                behavior) iteratively repacks and nudges per-type radii to match
                proportions. ``"graph_coloring"`` packs geometry once per
                replicate and assigns cell types via simulated-annealing graph
                coloring to match the target interaction statistics — far more
                consistent and faster-converging for interaction targets.
            coloring_params: Optional overrides for the ``"graph_coloring"`` SA
                schedule (keys: ``initial_temp``, ``final_temp``,
                ``cooling_rate``, ``max_iterations``, and optionally
                ``patience`` for adaptive stopping).
            n_restarts: For ``"graph_coloring"``, number of independent SA runs
                per replicate; the lowest-cost coloring is kept. >1 hardens
                against bad local minima at a linear cost. Default 1.
            radius_optimizer: For ``method="radius_tuning"``, the proportion
                tuner. ``"heuristic"`` (default, unchanged) uses the one-shot
                sqrt-ratio radius adjustment. ``"differential_evolution"`` runs
                a gradient-free scipy optimizer over per-type radius multipliers
                against a fixed-seed (deterministic) proportion objective.
                Gradient methods are deliberately not offered: the
                radius->cell-count map is integer-valued and stochastic, so
                finite-difference gradients are mostly zero. DE is more robust
                but slower; for matching interaction patterns prefer
                ``method="graph_coloring"``.
            de_params: Optional overrides for ``differential_evolution`` (e.g.
                ``maxiter``, ``popsize``, ``tol``).
        """
        if not NETWORKX_AVAILABLE:
            raise ImportError("NetworkX required for replicate generation")

        target_stats.validate()

        if method not in ("radius_tuning", "graph_coloring"):
            raise ValueError(
                f"method must be 'radius_tuning' or 'graph_coloring', got {method!r}."
            )
        if radius_optimizer not in ("heuristic", "differential_evolution"):
            raise ValueError(
                "radius_optimizer must be 'heuristic' or 'differential_evolution', "
                f"got {radius_optimizer!r}."
            )

        self.target_stats = target_stats
        self.tissue_dimensions = tissue_dimensions
        self.base_cell_radii = base_cell_radii
        self.network_mode = network_mode
        self.network_radius = network_radius
        self.seed = seed
        self.method = method
        self.n_restarts = max(1, int(n_restarts))
        self.radius_optimizer = radius_optimizer
        self.de_params = {'maxiter': 15, 'popsize': 10, 'tol': 0.01, 'polish': False}
        if de_params:
            self.de_params.update(de_params)

        # SA schedule for the graph-coloring path (overridable).
        self.coloring_params = {
            'initial_temp': 100.0,
            'final_temp': 0.1,
            'cooling_rate': 0.995,
            'max_iterations': 20000,
        }
        if coloring_params:
            self.coloring_params.update(coloring_params)

        # NOTE: we deliberately do NOT seed the global ``np.random`` module
        # RNG here. Doing so leaks state into the rest of the process and
        # makes reproducibility "best-effort" rather than guaranteed. Instead,
        # each replicate gets a deterministic per-replicate seed derived from
        # ``self.seed`` in ``generate_single_replicate``, and that seed is
        # threaded explicitly through ``TissueSection`` / ``SpherePacker``.

        # Extract cell types from target stats. Stored as a sorted tuple
        # rather than a set so iteration order is bit-stable across Python
        # processes (a plain set's order depends on PYTHONHASHSEED, and
        # downstream code feeds this iteration order into rng.choice via
        # dict construction — non-determinism there silently breaks the
        # reproducibility guarantee that ``seed`` is meant to provide).
        types_seen = set()
        for stat in target_stats.interaction_stats:
            types_seen.add(stat.type_a)
            types_seen.add(stat.type_b)
        self.cell_types = tuple(sorted(types_seen))

        # Set default proportions if not provided
        if target_stats.cell_type_proportions is None:
            n_types = len(self.cell_types)
            self.target_stats.cell_type_proportions = {
                ct: 1.0 / n_types for ct in self.cell_types
            }

        # Validate cell types match
        config_types = set(base_cell_radii.keys())
        missing = set(self.cell_types) - config_types
        if missing:
            raise ValueError(f"Cell types in target stats not in radii config: {missing}")
    
    def _compute_interaction_divergence(self,
                                       measured: List[InteractionStatistics],
                                       target: List[InteractionStatistics]) -> float:
        """
        Compute divergence between measured and target interaction statistics.

        Uses relative difference in normalized interactions as the metric.

        Per-pair semantics:
            - If both the target and measured value are zero, the pair has no
              signal in either direction; its contribution is ``nan`` rather
              than 0.0, so that "no signal" is not mistaken for a perfect
              match.
            - If the target value is > 0, the contribution is the relative
              difference ``|measured - target| / target``, capped at 2.0.
            - If the target value is 0 but the measured value is > 0, the
              contribution is the absolute difference (still capped at 2.0).
            - A pair present in ``target`` but missing from ``measured`` is
              treated as a 1.0 penalty (same as before).

        The aggregate score is the ``nanmean`` of the per-pair contributions,
        so pairs that were all-zero on both sides are ignored. If every pair
        is nan (the target has no signal at all), this returns ``nan``.

        Args:
            measured: Measured interaction statistics
            target: Target interaction statistics

        Returns:
            Average relative divergence (0 = perfect match, nan = no signal
            in any pair).
        """
        # Create lookup for measured stats
        measured_dict = {}
        for stat in measured:
            key = tuple(sorted([stat.type_a, stat.type_b]))
            measured_dict[key] = stat

        # Compute divergences
        divergences = []
        for target_stat in target:
            key = tuple(sorted([target_stat.type_a, target_stat.type_b]))

            if key not in measured_dict:
                # Missing interaction type - large penalty
                divergences.append(1.0)
                continue

            measured_stat = measured_dict[key]

            # Use normalized interactions as primary metric
            target_val = target_stat.normalized_interactions
            measured_val = measured_stat.normalized_interactions

            # All-zero pair: no signal anywhere, mark as nan so this pair
            # does NOT count as a perfect match in the aggregate.
            if target_val == 0 and measured_val == 0:
                divergences.append(np.nan)
                continue

            # Relative difference
            if target_val > 0:
                rel_diff = abs(measured_val - target_val) / target_val
            else:
                rel_diff = abs(measured_val - target_val)

            divergences.append(min(rel_diff, 2.0))  # Cap at 2.0

        if len(divergences) == 0:
            return float('nan')

        # nanmean ignores all-zero pairs; if every pair is nan (target has
        # no signal anywhere) numpy returns nan and emits a RuntimeWarning,
        # which we suppress because nan is the documented return value.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return float(np.nanmean(divergences))
    
    def _adjust_cell_type_proportions(self, 
                                     tissue: TissueSection,
                                     iteration: int) -> Dict[str, Tuple[float, float]]:
        """
        Adjust cell type proportions for next iteration.
        
        Args:
            tissue: Current tissue
            iteration: Current iteration number
        
        Returns:
            Updated cell radii configuration
        """
        # Get current proportions
        stats = tissue.get_cell_statistics()
        current_counts = stats['cell_types']
        total_cells = sum(current_counts.values())
        
        current_props = {
            ct: current_counts.get(ct, 0) / total_cells 
            for ct in self.cell_types
        }
        
        # Compute adjustment factors
        adjusted_radii = {}
        for cell_type in self.cell_types:
            target_prop = self.target_stats.cell_type_proportions[cell_type]
            current_prop = current_props[cell_type]
            
            # Adjust radii to change proportions
            # More cells -> smaller radii, fewer cells -> larger radii
            if current_prop > 0:
                adjustment = np.sqrt(target_prop / current_prop)
                adjustment = np.clip(adjustment, 0.7, 1.3)  # Limit adjustments
            else:
                adjustment = 1.1  # Slightly increase if missing
            
            min_r, max_r = self.base_cell_radii[cell_type]
            adjusted_radii[cell_type] = (
                min_r * adjustment,
                max_r * adjustment
            )

        return adjusted_radii

    @staticmethod
    def _round_proportions_to_counts(proportions: Dict[str, float],
                                     total: int) -> Dict[str, int]:
        """
        Round target proportions to integer counts summing exactly to ``total``.

        Uses the largest-remainder (Hamilton) method, with a stable
        (sorted-key) tie order so the result is reproducible across processes.
        """
        types = sorted(proportions.keys())
        if total <= 0:
            return {ct: 0 for ct in types}
        raw = {ct: proportions[ct] * total for ct in types}
        counts = {ct: int(np.floor(raw[ct])) for ct in types}
        remainder = total - sum(counts.values())
        # Hand the leftover units to the largest fractional parts.
        by_frac = sorted(types, key=lambda ct: (raw[ct] - counts[ct], ct), reverse=True)
        for ct in by_frac[:remainder]:
            counts[ct] += 1
        return counts

    def _build_colorizer_targets(self, graph) -> Dict:
        """
        Convert this generator's ``TargetStatistics`` into a GraphColorizer
        ``target_statistics`` dict for the *actual* packed ``graph``.

        Interaction targets are stored as ``normalized_interactions`` (fraction
        of possible edges). We de-normalize them against the packed graph's node
        count and the target cell-type proportions using the SAME
        ``max_possible`` convention as
        :meth:`SpatialNetworkAnalyzer.compute_interaction_statistics`
        (``n*(n-1)/2`` for self-pairs, ``na*nb`` for cross-pairs), so a perfect
        coloring of this geometry reproduces the target normalized interactions
        by construction.

        Returns:
            dict with ``node_counts`` (summing to the graph's node count),
            ``edge_counts`` (absolute counts, sorted ``'a-b'`` keys), and
            ``neighbor_dist`` (mean #c2-neighbors per c1-node).
        """
        N = graph.number_of_nodes()
        node_counts = self._round_proportions_to_counts(
            self.target_stats.cell_type_proportions, N
        )

        edge_counts: Dict[str, float] = {}
        neighbor_dist = defaultdict(lambda: defaultdict(float))
        for stat in self.target_stats.interaction_stats:
            a, b = sorted([stat.type_a, stat.type_b])
            na = node_counts.get(a, 0)
            nb = node_counts.get(b, 0)
            if a == b:
                max_possible = na * (na - 1) / 2 if na > 1 else 1
            else:
                max_possible = na * nb
            edges = stat.normalized_interactions * max_possible
            edge_counts['-'.join([a, b])] = edges

            # neighbor_dist[c1][c2] = expected #c2-neighbors per c1-node.
            # A same-type edge touches two c1 nodes, so it counts twice for that
            # type (matching GraphColorizer._calculate_statistics).
            if a == b:
                if na > 0:
                    neighbor_dist[a][a] = 2 * edges / na
            else:
                if na > 0:
                    neighbor_dist[a][b] = edges / na
                if nb > 0:
                    neighbor_dist[b][a] = edges / nb

        return {
            'node_counts': node_counts,
            'edge_counts': edge_counts,
            'neighbor_dist': neighbor_dist,
        }

    def _generate_single_replicate_colored(self,
                                           replicate_id: int,
                                           max_attempts: int = 1000,
                                           min_spacing: float = 0.5,
                                           allow_boundary: bool = True) -> Tuple[TissueSection, ReplicateStatistics]:
        """
        Generate one replicate via the graph-coloring method.

        Packs ONE tissue (fresh per-replicate seed -> geometric diversity), then
        assigns cell types with :func:`color_graph_to_targets` so the labeling
        matches the target interaction statistics on that fixed geometry. Unlike
        the radius-tuning path, the geometry is sampled once (not iterated); the
        simulated-annealing labeling is what gets optimized, which converges far
        more consistently.
        """
        # Per-replicate seed: identical convention to the radius-tuning path.
        if self.seed is not None:
            replicate_seed = int(
                np.random.SeedSequence([self.seed, replicate_id]).generate_state(1)[0]
            )
        else:
            replicate_seed = None

        # 1. Pack geometry once. Cell types from packing are placeholders; the
        #    coloring step overwrites them.
        tissue = TissueSection(
            height=self.tissue_dimensions[0],
            width=self.tissue_dimensions[1],
            thickness=self.tissue_dimensions[2],
            cell_radii=self.base_cell_radii,
            seed=replicate_seed,
        )
        num_cells = tissue.generate_cells(
            max_attempts=max_attempts,
            min_spacing=min_spacing,
            allow_boundary_cells=allow_boundary,
        )
        if num_cells == 0:
            raise RuntimeError(
                f"Replicate {replicate_id}: packing produced no cells"
            )

        # 2. Build the neighbor graph from the packed geometry.
        analyzer = SpatialNetworkAnalyzer()
        graph = analyzer.build_network_from_tissue(
            tissue, mode=self.network_mode, radius=self.network_radius
        )

        # 3. Derive GraphColorizer targets for THIS geometry, then color it.
        #    With n_restarts > 1, run independent SA attempts and keep the
        #    lowest-cost coloring (hardens against bad local minima).
        targets = self._build_colorizer_targets(graph)
        best_coloring = None
        best_cost = float('inf')
        for r in range(self.n_restarts):
            if self.n_restarts == 1:
                restart_seed = replicate_seed
            elif replicate_seed is not None:
                restart_seed = int(
                    np.random.SeedSequence([replicate_seed, r]).generate_state(1)[0]
                )
            else:
                restart_seed = None
            coloring, cost = color_graph_to_targets(
                graph,
                list(self.cell_types),
                targets,
                seed=restart_seed,
                return_cost=True,
                verbose=False,
                **self.coloring_params,
            )
            if cost < best_cost:
                best_cost = cost
                best_coloring = coloring
        coloring = best_coloring

        # 4. Write the coloring back onto the tissue cells and graph nodes so the
        #    measured statistics reflect the assigned types.
        for i, cell in enumerate(tissue.cells):
            if i in coloring:
                cell.cell_type = coloring[i]
        for node, color in coloring.items():
            graph.nodes[node]['cell_type'] = color

        # 5. Measure achieved interaction statistics and divergence.
        measured = analyzer.compute_interaction_statistics()
        divergence = self._compute_interaction_divergence(
            measured, self.target_stats.interaction_stats
        )

        tissue_stats = tissue.get_cell_statistics()
        replicate_stats = ReplicateStatistics(
            replicate_id=replicate_id,
            num_cells=tissue_stats['total_cells'],
            cell_type_counts=tissue_stats['cell_types'],
            packing_fraction=tissue_stats['packing_fraction'],
            interaction_stats=measured,
            divergence_score=divergence,
        )
        return tissue, replicate_stats

    def _generate_single_replicate_de(self,
                                      replicate_id: int,
                                      max_attempts: int = 1000,
                                      min_spacing: float = 0.5,
                                      allow_boundary: bool = True) -> Tuple[TissueSection, ReplicateStatistics]:
        """
        Radius-tuning replicate via gradient-free ``differential_evolution``.

        Optimizes per-type radius multipliers so the packed cell-type
        proportions match the target. The objective fixes the packing seed so
        it is deterministic for the optimizer (the radius->count map is
        otherwise stochastic and integer-valued, which defeats gradient
        methods). The solved radii are then used to pack and measure the
        replicate.
        """
        from scipy.optimize import differential_evolution

        if self.seed is not None:
            replicate_seed = int(
                np.random.SeedSequence([self.seed, replicate_id]).generate_state(1)[0]
            )
        else:
            replicate_seed = None

        types = list(self.cell_types)
        target_props = np.array(
            [self.target_stats.cell_type_proportions[t] for t in types]
        )

        def _radii_from_multipliers(mult):
            return {t: (self.base_cell_radii[t][0] * m, self.base_cell_radii[t][1] * m)
                    for t, m in zip(types, mult)}

        def _pack(radii):
            tissue = TissueSection(
                height=self.tissue_dimensions[0],
                width=self.tissue_dimensions[1],
                thickness=self.tissue_dimensions[2],
                cell_radii=radii,
                seed=replicate_seed,
            )
            tissue.generate_cells(
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary_cells=allow_boundary,
            )
            return tissue

        def objective(mult):
            tissue = _pack(_radii_from_multipliers(mult))
            stats = tissue.get_cell_statistics()
            total = stats['total_cells'] or 1
            props = np.array([stats['cell_types'].get(t, 0) / total for t in types])
            return float(np.sum((props - target_props) ** 2))

        bounds = [(0.5, 2.0)] * len(types)
        result = differential_evolution(
            objective, bounds, seed=replicate_seed, **self.de_params
        )

        # Final replicate with the solved radii.
        tissue = _pack(_radii_from_multipliers(result.x))
        analyzer = SpatialNetworkAnalyzer()
        analyzer.build_network_from_tissue(
            tissue, mode=self.network_mode, radius=self.network_radius
        )
        measured = analyzer.compute_interaction_statistics()
        divergence = self._compute_interaction_divergence(
            measured, self.target_stats.interaction_stats
        )

        tissue_stats = tissue.get_cell_statistics()
        replicate_stats = ReplicateStatistics(
            replicate_id=replicate_id,
            num_cells=tissue_stats['total_cells'],
            cell_type_counts=tissue_stats['cell_types'],
            packing_fraction=tissue_stats['packing_fraction'],
            interaction_stats=measured,
            divergence_score=divergence,
        )
        return tissue, replicate_stats

    def generate_single_replicate(self,
                                 replicate_id: int,
                                 max_attempts: int = 1000,
                                 min_spacing: float = 0.5,
                                 allow_boundary: bool = True,
                                 max_iterations: int = 5,
                                 tolerance: float = 0.15,
                                 patience: Optional[int] = None,
                                 method: Optional[str] = None) -> Tuple[TissueSection, ReplicateStatistics]:
        """
        Generate a single tissue replicate.

        For ``method="radius_tuning"`` (default), iteratively generates tissues
        and adjusts parameters to approach target statistics. For
        ``method="graph_coloring"``, packs geometry once and assigns cell types
        via simulated-annealing graph coloring (``max_iterations`` / ``tolerance``
        apply only to the radius-tuning loop).

        Args:
            replicate_id: Unique identifier for this replicate
            max_attempts: Max attempts for cell packing
            min_spacing: Minimum spacing between cells
            allow_boundary: Allow cells extending beyond bounds
            max_iterations: Max parameter adjustment iterations (radius-tuning only)
            tolerance: Acceptable divergence threshold (radius-tuning only)
            patience: Optional adaptive-stopping budget for the radius-tuning
                loop. When set, stop early once the best divergence has not
                improved for ``patience`` consecutive iterations. Lets you raise
                ``max_iterations`` without always paying for it.
            method: Override the instance ``method`` for this call.

        Returns:
            Tuple of (TissueSection, ReplicateStatistics)
        """
        if (method or self.method) == "graph_coloring":
            return self._generate_single_replicate_colored(
                replicate_id,
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary=allow_boundary,
            )

        if self.radius_optimizer == "differential_evolution":
            return self._generate_single_replicate_de(
                replicate_id,
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary=allow_boundary,
            )

        best_tissue = None
        best_divergence = float('inf')
        best_stats = None
        iters_since_improvement = 0

        current_radii = self.base_cell_radii.copy()

        # Derive a deterministic per-replicate seed from (self.seed, replicate_id)
        # using SeedSequence, which mixes the entropy in a stable, well-defined
        # way. When self.seed is None we fall through to None and tissue
        # generation remains unseeded (backward-compatible).
        if self.seed is not None:
            replicate_seed = int(
                np.random.SeedSequence([self.seed, replicate_id]).generate_state(1)[0]
            )
        else:
            replicate_seed = None

        for iteration in range(max_iterations):
            # Each iteration within a replicate gets its own derived seed so
            # the parameter-adjustment loop is also deterministic.
            if replicate_seed is not None:
                iter_seed = int(
                    np.random.SeedSequence([replicate_seed, iteration]).generate_state(1)[0]
                )
            else:
                iter_seed = None

            # Generate tissue
            tissue = TissueSection(
                height=self.tissue_dimensions[0],
                width=self.tissue_dimensions[1],
                thickness=self.tissue_dimensions[2],
                cell_radii=current_radii,
                seed=iter_seed,
            )

            num_cells = tissue.generate_cells(
                max_attempts=max_attempts,
                min_spacing=min_spacing,
                allow_boundary_cells=allow_boundary,
            )

            if num_cells == 0:
                warnings.warn(f"Replicate {replicate_id}, iteration {iteration}: No cells generated")
                continue

            # Analyze spatial interactions
            analyzer = SpatialNetworkAnalyzer()
            analyzer.build_network_from_tissue(
                tissue,
                mode=self.network_mode,
                radius=self.network_radius
            )

            measured_interactions = analyzer.compute_interaction_statistics()

            # Compute divergence
            divergence = self._compute_interaction_divergence(
                measured_interactions,
                self.target_stats.interaction_stats
            )

            # Track best result. nan divergence (no signal anywhere) is not
            # comparable; we still record it as the best if we have nothing.
            improved = False
            if best_tissue is None:
                best_divergence = divergence
                best_tissue = tissue
                best_stats = measured_interactions
            elif not np.isnan(divergence) and (np.isnan(best_divergence) or divergence < best_divergence):
                best_divergence = divergence
                best_tissue = tissue
                best_stats = measured_interactions
                improved = True

            # Check if we've met tolerance (nan never satisfies <= tolerance)
            if not np.isnan(divergence) and divergence <= tolerance:
                break

            # Adaptive stopping: bail out once the best divergence plateaus.
            iters_since_improvement = 0 if improved else iters_since_improvement + 1
            if patience is not None and iters_since_improvement >= patience:
                break

            # Adjust parameters for next iteration
            if iteration < max_iterations - 1:
                current_radii = self._adjust_cell_type_proportions(tissue, iteration)
        
        # Create statistics object
        if best_tissue is None:
            raise RuntimeError(f"Failed to generate replicate {replicate_id}")
        
        tissue_stats = best_tissue.get_cell_statistics()
        
        replicate_stats = ReplicateStatistics(
            replicate_id=replicate_id,
            num_cells=tissue_stats['total_cells'],
            cell_type_counts=tissue_stats['cell_types'],
            packing_fraction=tissue_stats['packing_fraction'],
            interaction_stats=best_stats,
            divergence_score=best_divergence
        )
        
        return best_tissue, replicate_stats
    
    def generate_replicates(self,
                           num_replicates: int,
                           max_attempts: int = 1000,
                           min_spacing: float = 0.5,
                           allow_boundary: bool = True,
                           max_iterations: int = 5,
                           tolerance: float = 0.15,
                           patience: Optional[int] = None,
                           parallel: bool = False,
                           max_workers: Optional[int] = None) -> List[Tuple[TissueSection, ReplicateStatistics]]:
        """
        Generate multiple tissue replicates.

        Args:
            num_replicates: Number of replicates to generate
            max_attempts: Max attempts for cell packing per replicate
            min_spacing: Minimum spacing between cells
            allow_boundary: Allow cells extending beyond bounds
            max_iterations: Max parameter adjustment iterations (radius-tuning only)
            tolerance: Acceptable divergence threshold (radius-tuning only)
            patience: Optional adaptive-stopping budget forwarded to each replicate.
            parallel: When True, generate replicates concurrently with a
                ``ProcessPoolExecutor``. Each replicate is independent and
                deterministically seeded from ``(self.seed, replicate_id)``, so
                results are identical to the serial path regardless of worker
                scheduling.
            max_workers: Worker count for the process pool (defaults to the
                executor's default when None).

        Returns:
            List of (TissueSection, ReplicateStatistics) tuples, in replicate order.
        """
        kwargs = dict(
            max_attempts=max_attempts,
            min_spacing=min_spacing,
            allow_boundary=allow_boundary,
            max_iterations=max_iterations,
            tolerance=tolerance,
            patience=patience,
        )

        if parallel and num_replicates > 1:
            from concurrent.futures import ProcessPoolExecutor, as_completed

            results: Dict[int, Tuple[TissueSection, ReplicateStatistics]] = {}
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(self.generate_single_replicate, replicate_id=i, **kwargs): i
                    for i in range(num_replicates)
                }
                for future in as_completed(futures):
                    i = futures[future]
                    results[i] = future.result()
                    stats = results[i][1]
                    print(f"Replicate {i + 1}/{num_replicates} done: "
                          f"Cells: {stats.num_cells}, Divergence: {stats.divergence_score:.4f}")
            return [results[i] for i in range(num_replicates)]

        replicates = []
        for i in range(num_replicates):
            print(f"Generating replicate {i+1}/{num_replicates}...")

            tissue, stats = self.generate_single_replicate(replicate_id=i, **kwargs)

            replicates.append((tissue, stats))
            print(f"  Cells: {stats.num_cells}, Divergence: {stats.divergence_score:.4f}")

        return replicates

    def consistency_report(self,
                           replicates,
                           metric: str = "divergence_score") -> Dict:
        """
        Quantify run-to-run consistency of generated replicates.

        Reports, per method, the mean / standard deviation / coefficient of
        variation of a chosen ReplicateStatistics metric, plus pairwise
        Cohen's d and the per-group N needed to distinguish methods (via
        :func:`power_analysis.compare_initialization_variance`). This is how
        you *measure* whether one strategy is more consistent than another.

        Note: a smaller coefficient of variation usually means more consistent,
        but CV is ``std/|mean|`` and inflates when the mean sits near zero
        (e.g. a method that matches the target almost perfectly). Read it
        alongside the absolute ``std`` and ``mean`` in the returned ``per_method``.

        Args:
            replicates: Either a list of ``(tissue, ReplicateStatistics)``
                tuples (a single method), or a dict mapping a method label to
                such a list (to compare methods, e.g.
                ``{"radius_tuning": [...], "graph_coloring": [...]}``).
            metric: ReplicateStatistics attribute to summarize
                (default ``"divergence_score"``).

        Returns:
            The dict from :func:`compare_initialization_variance`:
            ``per_method`` (n / mean / std / cv) and ``pairwise``
            (Cohen's d, required N per group). NaN metric values (pairs with no
            signal) are dropped before summarizing.
        """
        def _endpoints(reps):
            vals = [getattr(stats, metric) for _, stats in reps]
            return [v for v in vals if v is not None and not np.isnan(v)]

        if isinstance(replicates, dict):
            endpoints_by_method = {name: _endpoints(reps)
                                   for name, reps in replicates.items()}
        else:
            endpoints_by_method = {"replicates": _endpoints(replicates)}

        return compare_initialization_variance(endpoints_by_method)

    def export_replicate_statistics(self,
                                   replicates: List[Tuple[TissueSection, ReplicateStatistics]],
                                   base_filename: str):
        """
        Export statistics for all replicates to CSV files.
        
        Args:
            replicates: List of (tissue, stats) tuples
            base_filename: Base name for output files
        """
        # Summary statistics
        summary_data = []
        for tissue, stats in replicates:
            row = {
                'replicate_id': stats.replicate_id,
                'num_cells': stats.num_cells,
                'packing_fraction': stats.packing_fraction,
                'divergence_score': stats.divergence_score
            }
            # Add cell type counts
            for ct, count in stats.cell_type_counts.items():
                row[f'count_{ct}'] = count
                row[f'proportion_{ct}'] = count / stats.num_cells
            
            summary_data.append(row)
        
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv(f"{base_filename}_summary.csv", index=False)
        
        # Interaction statistics
        interaction_data = []
        for tissue, stats in replicates:
            for interaction in stats.interaction_stats:
                row = {
                    'replicate_id': stats.replicate_id,
                    'type_a': interaction.type_a,
                    'type_b': interaction.type_b,
                    'num_interactions': interaction.num_interactions,
                    'normalized_interactions': interaction.normalized_interactions,
                    'avg_distance': interaction.avg_distance,
                    'median_distance': interaction.median_distance
                }
                interaction_data.append(row)
        
        df_interactions = pd.DataFrame(interaction_data)
        df_interactions.to_csv(f"{base_filename}_interactions.csv", index=False)
        
        print(f"Exported statistics to {base_filename}_summary.csv and {base_filename}_interactions.csv")
    
    def export_replicate_tissues(self,
                                replicates: List[Tuple[TissueSection, ReplicateStatistics]],
                                output_dir: str):
        """
        Export each replicate tissue to a separate CSV file.
        
        Args:
            replicates: List of (tissue, stats) tuples
            output_dir: Directory for output files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for tissue, stats in replicates:
            filename = output_path / f"replicate_{stats.replicate_id:03d}_tissue.csv"
            tissue.export_to_csv(str(filename))
        
        print(f"Exported {len(replicates)} tissue files to {output_dir}")


def load_target_statistics_from_csv(filepath: str) -> TargetStatistics:
    """
    Load target statistics from a CSV file.
    
    Expected CSV format:
    type_a, type_b, num_interactions, normalized_interactions, avg_distance, median_distance
    
    Args:
        filepath: Path to CSV file
    
    Returns:
        TargetStatistics object
    """
    df = pd.read_csv(filepath)
    
    required_cols = ['type_a', 'type_b', 'normalized_interactions']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    
    # Create InteractionStatistics objects
    interaction_stats = []
    for _, row in df.iterrows():
        stat = InteractionStatistics(
            type_a=str(row['type_a']),
            type_b=str(row['type_b']),
            num_interactions=int(row.get('num_interactions', 0)),
            normalized_interactions=float(row['normalized_interactions']),
            avg_distance=float(row.get('avg_distance', 0.0)),
            median_distance=float(row.get('median_distance', 0.0))
        )
        interaction_stats.append(stat)
    
    return TargetStatistics(interaction_stats=interaction_stats)


def load_target_statistics_from_tissue(tissue: TissueSection,
                                      network_mode: str = "contact",
                                      network_radius: Optional[float] = None) -> TargetStatistics:
    """
    Extract target statistics from an existing tissue.
    
    Args:
        tissue: TissueSection to analyze
        network_mode: "contact" or "radius"
        network_radius: Distance threshold for "radius" mode
    
    Returns:
        TargetStatistics object
    """
    # Analyze the tissue
    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_tissue(
        tissue,
        mode=network_mode,
        radius=network_radius
    )
    
    # Get interaction statistics
    interaction_stats = analyzer.compute_interaction_statistics()
    
    # Get cell type proportions
    tissue_stats = tissue.get_cell_statistics()
    total_cells = tissue_stats['total_cells']
    cell_type_proportions = {
        ct: count / total_cells
        for ct, count in tissue_stats['cell_types'].items()
    }
    
    return TargetStatistics(
        interaction_stats=interaction_stats,
        cell_type_proportions=cell_type_proportions,
        target_cell_count=total_cells,
        target_density=tissue_stats['packing_fraction']
    )


def load_target_statistics_from_coordinates(filepath: str,
                                            network_mode: str = "contact",
                                            network_radius: Optional[float] = None) -> TargetStatistics:
    """
    Load FULL target statistics from a coordinate CSV file.

    This is a convenience composition: it reads a tissue from a coordinate
    CSV (one row per cell, with positions/radii/types) via
    ``load_tissue_from_csv`` and then derives target statistics from that
    reconstructed tissue via ``load_target_statistics_from_tissue``. It is
    exactly equivalent to::

        load_target_statistics_from_tissue(
            load_tissue_from_csv(filepath),
            network_mode=network_mode,
            network_radius=network_radius,
        )

    Because the statistics are computed from a real tissue, the returned
    ``TargetStatistics`` is fully populated: it includes the measured
    interaction statistics AND the ``cell_type_proportions`` and
    ``target_density`` (packing fraction) inferred from the cell coordinates.

    This is DISTINCT from ``load_target_statistics_from_csv``, which reads a
    precomputed interaction table (``type_a``, ``type_b``,
    ``normalized_interactions``, ...) and therefore does NOT populate
    ``cell_type_proportions`` or ``target_density``. Use this function when
    you have raw cell coordinates; use ``load_target_statistics_from_csv``
    when you already have an interaction-statistics table.

    Args:
        network_mode: "contact" (edges between touching cells) or "radius"
            (edges between cells within ``network_radius``), passed through
            to ``load_target_statistics_from_tissue``.
        network_radius: Distance threshold used when ``network_mode`` is
            "radius"; ignored otherwise.

    Returns:
        TargetStatistics object with interactions, cell type proportions,
        target cell count, and target density populated.
    """
    return load_target_statistics_from_tissue(
        load_tissue_from_csv(filepath),
        network_mode=network_mode,
        network_radius=network_radius,
    )
