"""Regression tests for `scripts/release.py` and the version strings it owns.

Two jobs:

1. **Guard the repo against version drift.** The version number lives in five
   files and only one of them is imported by anything, so the other four can
   disagree indefinitely without any runtime symptom. `test_repo_is_in_sync`
   turns that silent drift into a test failure.

2. **Guard the script's own edge cases.** The documented release flow writes
   the `## [X.Y.Z] - date` CHANGELOG section *before* running `bump`, so by
   then the newest heading already names the version being released. An
   earlier version of the script derived "previous version" from that heading,
   compared it to the target, found them equal, and silently skipped the
   compare link — the bug never showed up on the unhappy path, only on the
   documented one.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate
    raise RuntimeError("could not locate repository root")


@pytest.fixture(scope="module")
def release():
    """Load scripts/release.py as a module (it lives outside the package)."""
    path = _repo_root() / "scripts" / "release.py"
    spec = importlib.util.spec_from_file_location("release_script", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


CHANGELOG = """# Changelog

## [0.2.0] - 2026-08-01

### Added

- Something.

## [0.1.9] - 2026-07-01

### Fixed

- Something else.

[0.1.9]: https://example.invalid/compare/v0.1.8...v0.1.9
"""


@pytest.fixture
def fake_repo(tmp_path, release, monkeypatch):
    """Point the script at a throwaway CHANGELOG instead of the real repo."""
    (tmp_path / "CHANGELOG.md").write_text(CHANGELOG, encoding="utf-8")
    monkeypatch.setattr(release, "ROOT", tmp_path)
    return tmp_path


def test_repo_is_in_sync(release):
    """Every version string in the repo matches the newest CHANGELOG heading."""
    version, date = release.current_release()
    problems = release.audit(version, date)
    assert problems == [], (
        "Version strings have drifted from CHANGELOG.md:\n"
        + "\n".join(problems)
        + f"\n\nFix with: python scripts/release.py bump {version}"
    )


def test_previous_version_skips_the_version_being_released(fake_repo, release):
    """The regression: 0.2.0 is already the top heading when `bump 0.2.0` runs."""
    assert release.previous_version("0.2.0") == "0.1.9"


def test_previous_version_when_target_not_yet_in_changelog(fake_repo, release):
    assert release.previous_version("0.3.0") == "0.2.0"


def test_previous_version_is_none_for_a_first_release(tmp_path, release, monkeypatch):
    (tmp_path / "CHANGELOG.md").write_text("# Changelog\n", encoding="utf-8")
    monkeypatch.setattr(release, "ROOT", tmp_path)
    assert release.previous_version("0.1.0") is None


def test_compare_link_is_added_for_the_release_being_cut(fake_repo, release):
    release.add_compare_link("0.2.0", "0.1.9", dry_run=False)
    text = (fake_repo / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"[0.2.0]: {release.REPO_URL}/compare/v0.1.9...v0.2.0" in text
    # Newest first, and the existing reference survives.
    assert text.index("[0.2.0]: ") < text.index("[0.1.9]: ")


def test_compare_link_is_not_duplicated(fake_repo, release):
    release.add_compare_link("0.2.0", "0.1.9", dry_run=False)
    release.add_compare_link("0.2.0", "0.1.9", dry_run=False)
    text = (fake_repo / "CHANGELOG.md").read_text(encoding="utf-8")
    assert text.count("[0.2.0]: ") == 1


def test_dry_run_writes_nothing(fake_repo, release):
    before = (fake_repo / "CHANGELOG.md").read_text(encoding="utf-8")
    release.add_compare_link("0.2.0", "0.1.9", dry_run=True)
    assert (fake_repo / "CHANGELOG.md").read_text(encoding="utf-8") == before


def test_match_once_rejects_an_ambiguous_file(fake_repo, release):
    """A reformat that duplicates a version line must fail loudly, not silently."""
    import re

    (fake_repo / "dup.toml").write_text('version = "1"\nversion = "2"\n', encoding="utf-8")
    with pytest.raises(release.ReleaseError, match="expected exactly 1 match"):
        release._match_once("dup.toml", re.compile(r'^version = "([^"]+)"$', re.M), 'version = "1"\nversion = "2"\n')
