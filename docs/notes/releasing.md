# Release checklist

Maintainer-facing. Not part of the published docs site.

## The version-sync criterion

The version number appears in five files. Only one of them is imported by
anything, so the other four drift silently and nobody notices until a user
cites the wrong version or `pip show` disagrees with `__version__`:

| File | Form |
|---|---|
| `CHANGELOG.md` | `## [X.Y.Z] - YYYY-MM-DD` **(source of truth)** |
| `pyproject.toml` | `version = "X.Y.Z"` |
| `tissue_simulator/__init__.py` | `__version__ = "X.Y.Z"` |
| `CITATION.cff` | `version: X.Y.Z` and `date-released: "YYYY-MM-DD"` |
| `README.md` | BibTeX `version = {X.Y.Z}` and `year = {YYYY}` |

**Never edit these by hand.** The topmost `## [version] - date` heading in
`CHANGELOG.md` is authoritative; `scripts/release.py` syncs everything else to
it, so the release notes and the version bump can't disagree.

```bash
python scripts/release.py check          # verify sync — exit 1 if drifted
python scripts/release.py bump 0.1.16    # set the version everywhere
python scripts/release.py bump 0.1.16 --dry-run
```

`check` is cheap and side-effect free. Run it after any change that touches
version strings, packaging metadata, or citation metadata.

## Cutting a release

1. **Write the release notes.** Add a `## [X.Y.Z] - YYYY-MM-DD` section at the
   top of `CHANGELOG.md` describing what changed and *why* — the existing
   entries explain rationale, not just diffs. Match that depth.
2. **Bump.** `python scripts/release.py bump X.Y.Z` — from a clean working tree
   (the script refuses otherwise, so the bump lands as its own reviewable
   commit). This also appends the `[X.Y.Z]: …/compare/…` link reference.
3. **Test.** `pytest tests/ -v`
4. **Verify.** `python scripts/release.py check`
5. **Commit and open a PR.** Work lands on `main` through a pull request, as
   with #13/#14/#15:

   ```bash
   git commit -am "Release vX.Y.Z"
   git push -u origin <branch>
   gh pr create --base main
   ```

6. **Tag the merge commit on `main`,** after the PR merges:

   ```bash
   git checkout main && git pull
   git tag vX.Y.Z && git push --tags
   ```

   Tags are `vX.Y.Z`, matching every existing tag. Push the tag: Zenodo is
   wired to GitHub releases, and a missing tag breaks both the archive and the
   changelog compare link — `v0.1.7` was never tagged, which is why `0.1.8`'s
   link spans `v0.1.6...v0.1.8`. Note that `v0.1.15` and earlier point at
   *pre-merge branch commits*; tagging the merge commit on `main` is the
   intended convention going forward.
7. **Check Zenodo** minted a version DOI for the new tag. The *concept* DOI in
   `CITATION.cff` (`10.5281/zenodo.17465675`) always resolves to the latest
   release and does **not** change per release — leave it alone.

## Not on PyPI

`tissue_simulator` is **not published to PyPI**, so there is no build/upload
step. Users install from source. If that ever changes, add
`python -m build && twine upload dist/*` between the tag and the Zenodo check,
and note that a published version number can never be reused — a botched
upload burns that version.

## Citation metadata

`CITATION.cff` drives GitHub's "Cite this repository" button and the BibTeX
block in the README. If you change the author list, license, or DOI, change it
in `CITATION.cff` and `pyproject.toml` together — `release.py check` only
validates the version and date, not the other fields.
