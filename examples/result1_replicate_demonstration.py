"""
Result 1: Replicate tissues match a single target spatial-statistics fingerprint.

This headline demonstration shows that the ReplicateGenerator can produce
multiple distinct tissue samples whose spatial interaction statistics all
match a single reference fingerprint. We:

  1. Generate a deterministic reference tissue with three cell types
     (tumor, immune, stroma) in a 200x200x50 um box.
  2. Extract its spatial fingerprint with
     load_target_statistics_from_tissue (radius-mode network).
  3. Generate N_REPLICATES tissues with the ReplicateGenerator. Bit-for-bit
     reproducibility is now guaranteed for a fixed seed: TissueSection and
     SpherePacker accept an explicit seed argument that is threaded through
     to every random draw, and ReplicateGenerator derives a deterministic
     per-replicate seed from (seed, replicate_id) without touching the
     global numpy RNG.
  4. Measure each replicate's divergence vs. the target and write CSV
     summaries plus three matplotlib figures.

Run:
    python examples/result1_replicate_demonstration.py

Outputs are written to examples/output/result1/.
"""

import sys
import time
from pathlib import Path

# Make the package importable when running from a fresh checkout.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tissue_simulator import (
    TissueSection,
    ReplicateGenerator,
    SpatialNetworkAnalyzer,
    load_target_statistics_from_tissue,
)


# Default keeps the demo under 5 minutes on a laptop.
# The paper figure uses N_REPLICATES = 50.
N_REPLICATES = 10

REFERENCE_SEED = 20260513
TISSUE_DIMENSIONS = (200.0, 200.0, 50.0)  # (height, width, thickness) in um
CELL_RADII = {
    "tumor": (8.0, 12.0),
    "immune": (5.0, 8.0),
    "stroma": (6.0, 10.0),
}
NETWORK_MODE = "radius"
# Center-to-center distance threshold (um) for the spatial network. Chosen
# to be a small multiple of the typical cell radius so the contact graph
# is meaningfully populated regardless of MIN_SPACING.
NETWORK_RADIUS = 25.0
MAX_ATTEMPTS = 1500
MIN_SPACING = 0.5
ALLOW_BOUNDARY = True
MAX_ITERATIONS = 3
TOLERANCE = 0.15

OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "result1"

CELL_TYPE_COLORS = {
    "tumor": "#d62728",
    "immune": "#1f77b4",
    "stroma": "#2ca02c",
}


def generate_reference_tissue() -> TissueSection:
    """Generate a deterministic reference tissue."""
    tissue = TissueSection(
        height=TISSUE_DIMENSIONS[0],
        width=TISSUE_DIMENSIONS[1],
        thickness=TISSUE_DIMENSIONS[2],
        cell_radii=CELL_RADII,
        seed=REFERENCE_SEED,
    )
    n = tissue.generate_cells(
        max_attempts=MAX_ATTEMPTS,
        min_spacing=MIN_SPACING,
        allow_boundary_cells=ALLOW_BOUNDARY,
    )
    stats = tissue.get_cell_statistics()
    print(f"Reference tissue: {n} cells, "
          f"packing fraction = {stats['packing_fraction']:.4f}")
    print(f"  cell types: {stats['cell_types']}")
    return tissue


def edge_count_from_interactions(interactions) -> int:
    """Total number of edges across all pair-type interactions."""
    return int(sum(s.num_interactions for s in interactions))


def replicate_to_summary_row(rep_id: int, rep_stats, elapsed_s: float) -> dict:
    """Pack a ReplicateStatistics into a flat dict for CSV export."""
    row = {
        "replicate_id": rep_id,
        "n_cells": rep_stats.num_cells,
        "n_edges": edge_count_from_interactions(rep_stats.interaction_stats),
        "divergence_score": rep_stats.divergence_score,
        "packing_fraction": rep_stats.packing_fraction,
        "time_seconds": elapsed_s,
    }
    for ct, count in rep_stats.cell_type_counts.items():
        row[f"count_{ct}"] = count
    return row


def interactions_to_long_df(replicate_id, interactions) -> pd.DataFrame:
    """Convert a list of InteractionStatistics into a long-format DataFrame."""
    rows = []
    for s in interactions:
        pair = "-".join(sorted([s.type_a, s.type_b]))
        rows.append({
            "replicate_id": replicate_id,
            "pair": pair,
            "num_interactions": s.num_interactions,
            "normalized_interactions": s.normalized_interactions,
            "avg_distance": s.avg_distance,
            "median_distance": s.median_distance,
        })
    return pd.DataFrame(rows)


def plot_divergence_distribution(divergences, save_path: Path):
    """Histogram of replicate divergence scores with a smoothed density overlay."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    divergences = np.asarray(divergences, dtype=float)

    n_bins = max(5, min(15, len(divergences) // 2 + 3))
    counts, bin_edges, _ = ax.hist(
        divergences,
        bins=n_bins,
        color="#4c72b0",
        edgecolor="white",
        alpha=0.85,
        label="replicates",
    )

    # Lightweight KDE-like overlay using a Gaussian kernel.
    if len(divergences) >= 3 and divergences.std() > 0:
        xs = np.linspace(divergences.min(), divergences.max(), 200)
        bw = 1.06 * divergences.std() * len(divergences) ** (-1 / 5)
        bw = max(bw, 1e-4)
        kde = np.zeros_like(xs)
        for d in divergences:
            kde += np.exp(-0.5 * ((xs - d) / bw) ** 2)
        kde /= (len(divergences) * bw * np.sqrt(2 * np.pi))
        # Scale density to histogram counts for visual comparison.
        bin_width = bin_edges[1] - bin_edges[0]
        ax.plot(
            xs,
            kde * len(divergences) * bin_width,
            color="#c44e52",
            linewidth=2,
            label="kernel density",
        )

    ax.axvline(
        divergences.mean(),
        color="black",
        linestyle="--",
        linewidth=1,
        label=f"mean = {divergences.mean():.3f}",
    )
    ax.set_xlabel("Divergence vs. target")
    ax.set_ylabel("Replicate count")
    ax.set_title(f"Distribution of replicate divergence (N={len(divergences)})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)


def plot_replicate_examples(replicates, save_path: Path):
    """2x3 grid of xy-projection scatter plots for the first six replicates."""
    n_show = min(6, len(replicates))
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    seen_types = set()
    handles_by_type = {}

    for ax_idx in range(6):
        ax = axes_flat[ax_idx]
        if ax_idx >= n_show:
            ax.axis("off")
            continue
        tissue, rep_stats = replicates[ax_idx]
        xs, ys, colors, sizes = [], [], [], []
        for cell in tissue.cells:
            xs.append(cell.center[0])
            ys.append(cell.center[1])
            colors.append(CELL_TYPE_COLORS.get(cell.cell_type, "#888888"))
            # Marker area ~ radius^2 for a vaguely-correct visual scale.
            sizes.append((cell.radius ** 2) * 1.2)
        scatter = ax.scatter(
            xs, ys, c=colors, s=sizes, alpha=0.75, edgecolors="white", linewidths=0.3
        )
        ax.set_xlim(0, TISSUE_DIMENSIONS[1])
        ax.set_ylim(0, TISSUE_DIMENSIONS[0])
        ax.set_aspect("equal")
        ax.set_title(
            f"replicate {rep_stats.replicate_id}  "
            f"(n={rep_stats.num_cells}, div={rep_stats.divergence_score:.3f})",
            fontsize=10,
        )
        for cell in tissue.cells:
            if cell.cell_type not in seen_types:
                seen_types.add(cell.cell_type)
                handles_by_type[cell.cell_type] = plt.Line2D(
                    [0], [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=CELL_TYPE_COLORS.get(cell.cell_type, "#888888"),
                    markeredgecolor="white",
                    markersize=8,
                    label=cell.cell_type,
                )

    for ax in axes[-1, :]:
        ax.set_xlabel("x (um)")
    for ax in axes[:, 0]:
        ax.set_ylabel("y (um)")

    if handles_by_type:
        fig.legend(
            handles=list(handles_by_type.values()),
            loc="upper center",
            ncol=len(handles_by_type),
            bbox_to_anchor=(0.5, 1.02),
            frameon=False,
        )
    fig.suptitle("Replicate examples (xy projection of cell centers)", y=1.06)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_target_vs_achieved_bars(target_stats, replicates, save_path: Path):
    """Grouped bar chart comparing target vs. mean achieved edge counts per pair."""
    target_by_pair = {}
    for s in target_stats.interaction_stats:
        pair = "-".join(sorted([s.type_a, s.type_b]))
        target_by_pair[pair] = s.num_interactions

    achieved_by_pair = {pair: [] for pair in target_by_pair}
    for _, rep_stats in replicates:
        for s in rep_stats.interaction_stats:
            pair = "-".join(sorted([s.type_a, s.type_b]))
            if pair in achieved_by_pair:
                achieved_by_pair[pair].append(s.num_interactions)

    pairs = sorted(target_by_pair.keys())
    target_vals = [target_by_pair[p] for p in pairs]
    achieved_means = [np.mean(achieved_by_pair[p]) if achieved_by_pair[p] else 0.0
                      for p in pairs]
    achieved_stds = [np.std(achieved_by_pair[p]) if achieved_by_pair[p] else 0.0
                     for p in pairs]

    x = np.arange(len(pairs))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(7, 1.2 * len(pairs) + 3), 5))
    ax.bar(x - width / 2, target_vals, width,
           label="target", color="#2c3e50", alpha=0.9)
    ax.bar(x + width / 2, achieved_means, width,
           yerr=achieved_stds, capsize=4,
           label=f"replicate mean (N={len(replicates)})",
           color="#e67e22", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, rotation=30, ha="right")
    ax.set_ylabel(f"Edge count ({NETWORK_MODE}-mode interactions)")
    ax.set_title("Target vs. achieved edge counts per cell-type pair")
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    overall_start = time.time()

    # 1. Reference tissue.
    print("[1/5] Generating reference tissue...")
    reference_tissue = generate_reference_tissue()

    # 2. Target fingerprint.
    print("[2/5] Extracting target statistics from reference tissue...")
    target_stats = load_target_statistics_from_tissue(
        reference_tissue,
        network_mode=NETWORK_MODE,
        network_radius=NETWORK_RADIUS,
    )
    # Also pull reference edge count via SpatialNetworkAnalyzer so we can
    # report it alongside the per-replicate edge counts.
    ref_analyzer = SpatialNetworkAnalyzer()
    ref_analyzer.build_network_from_tissue(
        reference_tissue, mode=NETWORK_MODE, radius=NETWORK_RADIUS
    )
    ref_global = ref_analyzer.compute_global_statistics()
    print(f"  target interaction types: {len(target_stats.interaction_stats)}")
    print(f"  target n_cells = {target_stats.target_cell_count}, "
          f"target n_edges = {ref_global.total_edges}")
    if ref_global.total_edges == 0:
        raise RuntimeError(
            "Reference network has zero edges; the divergence comparison "
            "would be vacuous. Increase NETWORK_RADIUS, switch to "
            "network_mode='contact' with a small MIN_SPACING, or expand "
            "the tissue dimensions so cells form an interaction network."
        )

    # 3. Generate replicates.
    print(f"[3/5] Generating {N_REPLICATES} replicates...")
    replicates = []
    summary_rows = []
    achieved_long = []

    for i in range(N_REPLICATES):
        # Each replicate uses a deterministic per-run seed. The generator
        # threads this seed through TissueSection -> SpherePacker without
        # touching the global numpy RNG, so the run is bit-for-bit
        # reproducible.
        generator = ReplicateGenerator(
            target_stats=target_stats,
            tissue_dimensions=TISSUE_DIMENSIONS,
            base_cell_radii=CELL_RADII,
            network_mode=NETWORK_MODE,
            network_radius=NETWORK_RADIUS,
            seed=REFERENCE_SEED + 1 + i,
        )
        t0 = time.time()
        tissue, rep_stats = generator.generate_single_replicate(
            replicate_id=i,
            max_attempts=MAX_ATTEMPTS,
            min_spacing=MIN_SPACING,
            allow_boundary=ALLOW_BOUNDARY,
            max_iterations=MAX_ITERATIONS,
            tolerance=TOLERANCE,
        )
        elapsed = time.time() - t0

        replicates.append((tissue, rep_stats))
        summary_rows.append(replicate_to_summary_row(i, rep_stats, elapsed))
        achieved_long.append(interactions_to_long_df(i, rep_stats.interaction_stats))

        print(f"  replicate {i}: n_cells={rep_stats.num_cells}, "
              f"divergence={rep_stats.divergence_score:.4f}, "
              f"time={elapsed:.1f}s")

    # 4. Write CSV outputs.
    print("[4/5] Writing CSV outputs...")
    summary_df = pd.DataFrame(summary_rows)
    summary_path = OUTPUT_DIR / "replicate_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    # Target vs achieved interaction-level table (one row per pair per replicate
    # plus the target row, replicate_id = -1).
    target_rows = []
    for s in target_stats.interaction_stats:
        pair = "-".join(sorted([s.type_a, s.type_b]))
        target_rows.append({
            "replicate_id": -1,  # -1 denotes the reference target
            "pair": pair,
            "num_interactions": s.num_interactions,
            "normalized_interactions": s.normalized_interactions,
            "avg_distance": s.avg_distance,
            "median_distance": s.median_distance,
        })
    achieved_df = pd.concat(
        [pd.DataFrame(target_rows)] + achieved_long,
        ignore_index=True,
    )
    achieved_path = OUTPUT_DIR / "target_vs_achieved.csv"
    achieved_df.to_csv(achieved_path, index=False)

    # 5. Figures.
    print("[5/5] Plotting figures...")
    divergences = [r.divergence_score for _, r in replicates]
    plot_divergence_distribution(divergences,
                                 OUTPUT_DIR / "divergence_distribution.png")
    plot_replicate_examples(replicates,
                            OUTPUT_DIR / "replicate_examples.png")
    plot_target_vs_achieved_bars(target_stats, replicates,
                                 OUTPUT_DIR / "target_vs_achieved_bar.png")

    # Concise summary.
    total_elapsed = time.time() - overall_start
    n_cells_array = np.array([r.num_cells for _, r in replicates])
    div_array = np.array(divergences, dtype=float)
    n_nan_div = int(np.isnan(div_array).sum())
    print()
    print("=" * 60)
    print("Result 1 summary")
    print("=" * 60)
    print(f"  Target n_cells:          {target_stats.target_cell_count}")
    print(f"  Achieved n_cells:        "
          f"mean={n_cells_array.mean():.1f}, std={n_cells_array.std():.1f}")
    print(f"  Mean divergence:         "
          f"{np.nanmean(div_array):.4f} (std {np.nanstd(div_array):.4f})")
    if n_nan_div > 0:
        print(f"  WARNING: {n_nan_div} replicate(s) had nan divergence "
              f"(target had no signal in any pair); excluded from the mean.")
    print(f"  Replicates generated:    {len(replicates)}")
    print(f"  Total run time:          {total_elapsed:.1f} s")
    print(f"  Output directory:        {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
