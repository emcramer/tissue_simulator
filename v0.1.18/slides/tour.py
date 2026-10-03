import marimo

__generated_with = "0.19.11"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    return (mo,)


@app.cell
def __(mo):
    mo.md(
        r"""
        # tissue_simulator: an end-to-end tour

        3D simulated biological tissue sections with network-based spatial analysis
        and ABM-initialization workflows.

        [Docs](https://emcramer.github.io/tissue_simulator/) · [Source](https://github.com/emcramer/tissue_simulator) · [Wiki](https://github.com/emcramer/tissue_simulator/wiki)
        """
    )
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        ## What it does

        - **Pack** 3D tissues with random sphere placement (`TissueSection` + `SpherePacker`).
        - **Slice** them at arbitrary planes (`TissueSlicer`).
        - **Build spatial networks** and compute interaction statistics (`SpatialNetworkAnalyzer`).
        - **Generate replicates** matching a target's contact statistics (`ReplicateGenerator`).
        - **Assign cell types** by simulated annealing to match an adjacency target (`GraphColorizer`).
        - **Reproducibility** — every stochastic component takes a `seed=`.
        - **MCP server** exposes the API to LLMs (Claude Desktop / agents).
        """
    )
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Install

        ```bash
        pip install -e ".[mcp]"
        ```

        Optional extras: `.[docs]` for the docs site, `.[dev]` for tests.
        """
    )
    return


@app.cell
def __():
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt

    # Marimo renders on a scrolling page at natural size, so modest defaults
    # give nice aspect ratios without needing per-cell resizing.
    plt.rcParams["figure.figsize"] = (8, 5)
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["savefig.dpi"] = 100
    plt.rcParams["savefig.bbox"] = "tight"

    import tissue_simulator
    print("tissue_simulator", tissue_simulator.__version__)
    return (plt,)


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Create a tissue

        A few lines packs a 3D tissue. `seed=2026` makes it reproducible.
        """
    )
    return


@app.cell
def __(plt):
    from tissue_simulator import TissueSection

    tissue = TissueSection(height=150, width=150, thickness=40, cell_radii=(5, 9), seed=2026)
    n = tissue.generate_cells(max_attempts=500)
    print(f"Placed {n} cells; packing fraction = {tissue.get_cell_statistics()['packing_fraction']:.2f}")
    tissue.visualize(elevation=22, azimuth=35)
    plt.gca()
    return TissueSection, n, tissue


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Multi-type tissues

        Pass a dict of `(min, max)` radii per type. Cells are color-coded.
        """
    )
    return


@app.cell
def __(TissueSection, plt):
    multi = TissueSection(
        height=150, width=150, thickness=40,
        cell_radii={"cancer": (5, 8), "immune": (3, 6), "stroma": (8, 12)},
        seed=2026,
    )
    multi.generate_cells(max_attempts=500)
    multi.visualize(elevation=22, azimuth=35)
    plt.gca()
    return (multi,)


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Slice it (2D)

        `TissueSlicer` extracts a planar cross-section. Cell radii project to circles.
        """
    )
    return


@app.cell
def __(multi, plt):
    from tissue_simulator import TissueSlicer

    slicer = TissueSlicer(multi)
    slicer.slice_plane(z_position=20)
    slicer.visualize_slice_2d(
        title=f"z = 20 μm slice ({len(slicer.slice_cells)} cells)",
        figsize=(8, 5),
    )
    plt.gca()
    return TissueSlicer, slicer


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Slice in 3D context

        Same slice, shown against the parent tissue.
        """
    )
    return


@app.cell
def __(plt, slicer):
    slicer.visualize_slice_in_3d(plane_alpha=0.25)
    plt.gca()
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Build a spatial network

        Connect cells whose centers are within a radius — here 18 μm.
        """
    )
    return


@app.cell
def __(multi, plt):
    from tissue_simulator import SpatialNetworkAnalyzer

    analyzer = SpatialNetworkAnalyzer()
    analyzer.build_network_from_tissue(multi, mode="radius", radius=18.0)
    analyzer.visualize_network(layout="spatial", figsize=(8, 5))
    plt.gca()
    return SpatialNetworkAnalyzer, analyzer


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Per-type interaction statistics

        How many edges of each cross-type pair? Normalised by the maximum possible.
        """
    )
    return


@app.cell
def __(analyzer, plt):
    stats = analyzer.compute_interaction_statistics()
    labels = [f"{s.type_a}–{s.type_b}" for s in stats]
    vals = [s.normalized_interactions for s in stats]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(labels, vals)
    ax.set_xlabel("normalized interactions")
    ax.set_title("Cross-type contact frequency")
    fig.tight_layout()
    fig
    return ax, fig, labels, stats, vals


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Generate replicates

        `ReplicateGenerator` repacks new positions that match the target's contact statistics
        and proportions. Same `seed=` → same replicates.
        """
    )
    return


@app.cell
def __(multi, plt):
    from tissue_simulator import ReplicateGenerator, load_target_statistics_from_tissue

    target_stats = load_target_statistics_from_tissue(multi, network_mode="radius", network_radius=18.0)
    gen = ReplicateGenerator(
        target_stats=target_stats,
        tissue_dimensions=(150, 150, 40),
        base_cell_radii={"cancer": (5, 8), "immune": (3, 6), "stroma": (8, 12)},
        network_mode="radius",
        network_radius=18.0,
        seed=42,
    )
    replicates = gen.generate_replicates(num_replicates=2, max_attempts=300, max_iterations=2, tolerance=0.25)
    print(
        f"Generated {len(replicates)} replicates, "
        f"{len(replicates[0][0].cells)} and {len(replicates[1][0].cells)} cells respectively."
    )
    replicates[0][0].visualize(elevation=22, azimuth=35)
    plt.gca()
    return (
        ReplicateGenerator,
        gen,
        load_target_statistics_from_tissue,
        replicates,
        target_stats,
    )


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Cell-type assignment via simulated annealing

        `GraphColorizer` matches a target adjacency structure on a fixed graph. Seeded for reproducibility.
        """
    )
    return


@app.cell
def __(plt):
    import networkx as nx
    from tissue_simulator import GraphColorizer
    from tissue_simulator.graph_coloring import visualize_colored_graph

    # Small synthetic graph for the demo — quick and visually clean.
    G = nx.connected_watts_strogatz_graph(40, 4, 0.2, seed=7)
    target = {
        "node_counts": {"cancer": 14, "immune": 14, "stroma": 12},
        "edge_counts": {
            "cancer-cancer": 12, "cancer-immune": 14, "cancer-stroma": 10,
            "immune-immune": 12, "immune-stroma": 10, "stroma-stroma": 8,
        },
        "neighbor_dist": {
            "cancer": {"cancer": 1.5, "immune": 1.3, "stroma": 1.0},
            "immune": {"cancer": 1.3, "immune": 1.5, "stroma": 1.0},
            "stroma": {"cancer": 1.0, "immune": 1.0, "stroma": 1.4},
        },
    }
    colorizer = GraphColorizer(
        target_graph=G, colors=["cancer", "immune", "stroma"],
        target_statistics=target, seed=42,
    )
    assignment = colorizer.colorize(
        initial_temp=10.0, final_temp=0.5,
        cooling_rate=0.95, max_iterations=200, verbose=False,
    )
    visualize_colored_graph(G, assignment, title="GraphColorizer result", figsize=(8, 5))
    plt.gca()
    return (
        G,
        GraphColorizer,
        assignment,
        colorizer,
        nx,
        target,
        visualize_colored_graph,
    )


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Reproducibility (`seed=`)

        Two runs with the same seed produce identical outputs.
        """
    )
    return


@app.cell
def __(G, GraphColorizer, target):
    def run_once(seed):
        g = GraphColorizer(
            target_graph=G, colors=["cancer", "immune", "stroma"],
            target_statistics=target, seed=seed,
        )
        return g.colorize(
            initial_temp=10.0, final_temp=0.5, cooling_rate=0.95,
            max_iterations=200, verbose=False,
        )

    a1, a2 = run_once(42), run_once(42)
    b = run_once(43)
    print(f"Same seed equal? {a1 == a2}")
    print(f"Different seed equal? {a1 == b}")
    return a1, a2, b, run_once


@app.cell
def __(mo):
    mo.md(
        r"""
        ## Where next

        - 📚 **Docs site:** <https://emcramer.github.io/tissue_simulator/>
        - 🛠️ **Source / issues:** <https://github.com/emcramer/tissue_simulator>
        - 🤖 **MCP guide:** <https://emcramer.github.io/tissue_simulator/latest/guides/mcp/>
        - 💬 **Wiki (FAQ / Troubleshooting / Roadmap):** <https://github.com/emcramer/tissue_simulator/wiki>

        Thanks — bug reports + feature ideas welcome on the issue tracker.
        """
    )
    return


if __name__ == "__main__":
    app.run()
