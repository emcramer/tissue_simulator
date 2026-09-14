"""
Density-aware replicates on a synthetic tumor-nest process.

A reference process with dense tumor nests, an immune margin and sparse stroma
is sampled independently of the replicate machinery. For a few target samples,
replicates are built on the uniform scaffold (cell count calibrated to the
target) and on the density-aware scaffold, and both are compared with the
natural variability of the process.

Held-out statistics, matched by neither method: 9 k-nearest-neighbor
composition entries and 36 cross-type L(r) values (6 type pairs at
r = 10-80 µm, border corrected). |z| is the distance of the replicate mean
from the truth mean in truth standard deviations. The "target itself" row is
the floor: replicates of one sample cannot be closer to the process than the
sample is.

Run:
    python examples/density_aware_replicates.py          # 20 truth samples, 3 targets x 10 replicates
    python examples/density_aware_replicates.py --quick  # 8 truth samples, 1 target x 4 replicates
"""

import argparse
import contextlib
import io
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).parent.parent))

from tissue_simulator import DensityModel, ReplicateGenerator, SpherePacker, TissueSection
from tissue_simulator.replicate_generator import load_target_statistics_from_tissue

SIZE = 400.0
THICKNESS = 1.0
TYPES = ("CD8", "Stroma", "Tumor")
L_RADII = (10, 20, 30, 40, 60, 80)
RADII = {t: (3.5, 5.0) for t in TYPES}
COLORING = dict(cooling_rate=0.9995, max_iterations=20000)
GROUPS = {
    "kNN": list(range(9)),
    "L r<=20": [9 + 6 * p + i for p in range(6) for i in (0, 1)],
    "L 30-40": [9 + 6 * p + i for p in range(6) for i in (2, 3)],
    "L 60-80": [9 + 6 * p + i for p in range(6) for i in (4, 5)],
}


def nest_process(seed):
    """Random disc nests at full density, a 12 µm margin at 70%, stroma at 30%."""
    rng = np.random.default_rng(seed)
    packed = SpherePacker((SIZE, SIZE, THICKNESS), {"c": (3.5, 5.0)}, min_spacing=0.3,
                          seed=int(rng.integers(2**31))).pack(max_attempts=1000)
    centers = rng.uniform(0, SIZE, (5, 2))
    radii = rng.uniform(35, 55, 5)
    tissue = TissueSection(SIZE, SIZE, THICKNESS, RADII)
    for cell in packed:
        d = (np.linalg.norm(centers - cell.center[:2], axis=1) - radii).min()
        if rng.random() < (1.0 if d < 0 else 0.7 if d < 12 else 0.3):
            u = rng.random()
            if d < 0:
                cell.cell_type = "Tumor" if u < 0.9 else "CD8"
            elif d < 12:
                cell.cell_type = "CD8" if u < 0.7 else "Stroma"
            else:
                cell.cell_type = "Stroma" if u < 0.8 else "CD8"
            tissue.cells.append(cell)
    return tissue


def _xy_types(tissue):
    return (np.array([c.center[:2] for c in tissue.cells]),
            np.array([c.cell_type for c in tissue.cells]))


def held_out_statistics(tissue):
    """9 kNN-composition entries followed by 36 border-corrected cross-L values."""
    xy, types = _xy_types(tissue)
    tree = cKDTree(xy)
    _, neighbors = tree.query(xy, k=11)
    stats = []
    for a in TYPES:
        sel = types == a
        stats += [float((types[neighbors[sel, 1:]] == b).mean()) if sel.any() else 0.0
                  for b in TYPES]
    edge = np.minimum.reduce([xy[:, 0], xy[:, 1], SIZE - xy[:, 0], SIZE - xy[:, 1]])
    for i, a in enumerate(TYPES):
        for b in TYPES[i:]:
            focal_all = np.flatnonzero(types == a)
            is_b = types == b
            n_b = is_b.sum() - (a == b)
            for r in L_RADII:
                focal = focal_all[edge[focal_all] > r]
                if focal.size == 0 or n_b <= 0:
                    stats.append(0.0)
                    continue
                count = sum(int(is_b[ball].sum()) - (a == b)
                            for ball in tree.query_ball_point(xy[focal], r))
                stats.append(float(np.sqrt(SIZE * SIZE * count / (focal.size * n_b) / np.pi)))
    return np.array(stats)


def mean_degree(tissue, radius=20.0):
    xy, _ = _xy_types(tissue)
    return 2 * len(cKDTree(xy).query_pairs(radius)) / len(xy)


def shared_positions(a, b, tolerance=3.0):
    """Fraction of the smaller configuration matched one-to-one within ``tolerance``."""
    xa, _ = _xy_types(a)
    xb, _ = _xy_types(b)
    if len(xa) > len(xb):
        xa, xb = xb, xa
    cost = np.linalg.norm(xa[:, None] - xb[None], axis=2)
    rows, cols = linear_sum_assignment(np.minimum(cost, 10 * tolerance))
    return float((cost[rows, cols] <= tolerance).mean())


def calibrated_max_attempts(n_target, seed):
    """Uniform-scaffold threshold whose packed count is closest to ``n_target``."""
    candidates = [2, 3, 5, 8, 12, 20, 35, 60, 100, 200, 400]
    counts = [len(SpherePacker((SIZE, SIZE, THICKNESS), RADII, seed=seed).pack(max_attempts=m))
              for m in candidates]
    return candidates[int(np.argmin(np.abs(np.array(counts) - n_target)))]


def build_replicates(target_tissue, scaffold, n, seed, composition_bin=40.0):
    target = load_target_statistics_from_tissue(target_tissue, network_mode="radius",
                                                network_radius=20.0)
    target.target_density = None  # not meaningful for a 1 µm slab
    kwargs = {}
    max_attempts = 1000
    if scaffold == "density-aware":
        kwargs["density_model"] = DensityModel.from_tissue(target_tissue, seed=seed)
        kwargs["composition_bin"] = composition_bin
    else:
        max_attempts = calibrated_max_attempts(len(target_tissue.cells), seed)
    generator = ReplicateGenerator(target, (SIZE, SIZE, THICKNESS), RADII,
                                   network_mode="radius", network_radius=20.0, seed=seed,
                                   method="graph_coloring", coloring_params=COLORING, **kwargs)
    tissues, start = [], time.time()
    for i in range(n):
        with contextlib.redirect_stdout(io.StringIO()):
            tissues.append(generator.generate_single_replicate(i, max_attempts=max_attempts)[0])
    return tissues, (time.time() - start) / n


def median_pairwise_distance(matrix, scale):
    z = matrix / scale
    return float(np.median([np.linalg.norm(z[i] - z[j])
                            for i in range(len(z)) for j in range(i + 1, len(z))]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--quick", action="store_true", help="fewer samples, for a smoke run")
    parser.add_argument("--composition-bin", type=float, default=40.0,
                        help="side of the annealer's composition bins in µm")
    args = parser.parse_args()
    n_truth, n_reps, n_targets = (8, 4, 1) if args.quick else (20, 10, 3)

    truth_tissues = [nest_process(s) for s in range(n_truth)]
    truth = np.array([held_out_statistics(t) for t in truth_tissues])
    mu, sd = truth.mean(axis=0), truth.std(axis=0, ddof=1)
    sd[sd < 1e-9] = 1e-9
    truth_spread = median_pairwise_distance(truth, sd)
    chance = np.mean([shared_positions(truth_tissues[i], truth_tissues[i + 1])
                      for i in range(min(4, n_truth - 1))])

    rows = {"target itself": [], "uniform": [], "density-aware": []}
    for t in range(n_targets):
        target = truth_tissues[t]
        rows["target itself"].append({"z": np.abs((truth[t] - mu) / sd)})
        for scaffold in ("uniform", "density-aware"):
            tissues, seconds = build_replicates(target, scaffold, n_reps, seed=1 + t,
                                                composition_bin=args.composition_bin)
            reps = np.array([held_out_statistics(x) for x in tissues])
            rows[scaffold].append({
                "z": np.abs((reps.mean(axis=0) - mu) / sd),
                "spread": median_pairwise_distance(reps, sd) / truth_spread,
                "shared": np.mean([shared_positions(tissues[i], tissues[i + 1])
                                   for i in range(n_reps - 1)]),
                "degree": np.mean([mean_degree(x) for x in tissues]) / mean_degree(target),
                "seconds": seconds,
            })

    print(f"{n_truth} truth samples, {n_targets} target(s) x {n_reps} replicates; "
          f"chance shared positions {chance:.0%}\n")
    header = f"{'':<14} {'mean|z|':>7} {'|z|>1':>7} " + " ".join(f"{g:>8}" for g in GROUPS)
    print(header + f" {'spread':>7} {'shared':>7} {'degree':>7} {'s/rep':>6}")
    for name, entries in rows.items():
        z = np.array([e["z"] for e in entries])
        line = (f"{name:<14} {z.mean():7.2f} {np.mean((z > 1).sum(axis=1)):5.1f}/45 "
                + " ".join(f"{z[:, idx].mean():8.2f}" for idx in GROUPS.values()))
        if "spread" in entries[0]:
            line += " ".join(["", *(f"{np.mean([e[k] for e in entries]):7.2f}"
                                    for k in ("spread", "shared", "degree")),
                              f"{np.mean([e['seconds'] for e in entries]):6.1f}"])
        print(line)


if __name__ == "__main__":
    main()
