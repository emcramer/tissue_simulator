# scripts/

Helper scripts kept out of the importable package so they don't pollute
`pip install tissue_simulator`.

| Script | What it does |
|---|---|
| `configure_mcp.sh` | Writes a Claude Desktop config block that points at this checkout's `run_mcp_server.py`. Idempotent; backs up any existing config. |
| `release.py` | Keeps the version string in sync across `CHANGELOG.md`, `pyproject.toml`, `__init__.py`, `CITATION.cff`, and the README BibTeX. `check` verifies, `bump X.Y.Z` sets. See `docs/notes/releasing.md`. |

Run them from the repository root, e.g. `./scripts/configure_mcp.sh`.
