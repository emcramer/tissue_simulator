<!-- mermaid-ai-skills:start -->
## Mermaid Diagrams

When the user asks to create, edit, or visualize a diagram, follow the
instructions in `.github/instructions/mermaid.instructions.md`.
<!-- mermaid-ai-skills:end -->

## Versioning

The version string lives in five files: `CHANGELOG.md` (source of truth),
`pyproject.toml`, `tissue_simulator/__init__.py`, `CITATION.cff`, and the
BibTeX block in `README.md`.

Never hand-edit a version string — they drift. Use `scripts/release.py`:

```bash
python scripts/release.py check        # verify sync; exit 1 if drifted
python scripts/release.py bump 0.1.16  # set it everywhere
```

Run `check` after touching packaging or citation metadata. Do not bump the
version as a side effect of a feature change. Full procedure:
`docs/notes/releasing.md`.
