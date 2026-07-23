#!/usr/bin/env python3
"""Keep every version string in the repo in sync, and drive a release.

The version number lives in five places. They drift silently, because nothing
imports four of them:

    pyproject.toml            version = "X.Y.Z"          <- what pip installs
    tissue_simulator/__init__ __version__ = "X.Y.Z"      <- what users see
    CITATION.cff              version / date-released    <- what GitHub cites
    README.md                 BibTeX version / year      <- what people paste
    CHANGELOG.md              ## [X.Y.Z] - YYYY-MM-DD    <- source of truth

`CHANGELOG.md` is authoritative: its topmost `## [version] - date` heading is
what every other site is synced *to*. That keeps the release notes and the
version bump from ever disagreeing, and means this script never has to invent
a version or a date.

Usage
-----
    python scripts/release.py check              # verify sync; exit 1 if drifted
    python scripts/release.py bump 0.1.16        # bump everything to 0.1.16
    python scripts/release.py bump 0.1.16 --date 2026-07-23
    python scripts/release.py bump 0.1.16 --dry-run

`check` is safe to run any time and is the thing CI (or an agent finishing a
change) should run. `bump` rewrites files; run it from a clean working tree.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO_URL = "https://github.com/emcramer/tissue_simulator"

SEMVER = re.compile(r"^\d+\.\d+\.\d+([-+.][0-9A-Za-z.\-+]+)?$")
ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

# Latest release heading in CHANGELOG.md — the source of truth.
CHANGELOG_HEADING = re.compile(r"^## \[(?P<v>[^\]]+)\] - (?P<d>\d{4}-\d{2}-\d{2})$", re.M)

# (relative path, regex with a named group, replacement template).
# Each pattern MUST match exactly once; a zero- or multi-match is a hard error
# so that a reformat upstream fails loudly instead of silently skipping a file.
VERSION_SITES = [
    ("pyproject.toml", re.compile(r'^version = "(?P<v>[^"]+)"$', re.M), 'version = "{v}"'),
    (
        "tissue_simulator/__init__.py",
        re.compile(r'^__version__ = "(?P<v>[^"]+)"$', re.M),
        '__version__ = "{v}"',
    ),
    ("CITATION.cff", re.compile(r"^version: (?P<v>\S+)$", re.M), "version: {v}"),
    ("README.md", re.compile(r"^  version   = \{(?P<v>[^}]+)\},$", re.M), "  version   = {{{v}}},"),
]

DATE_SITES = [
    (
        "CITATION.cff",
        re.compile(r'^date-released: "(?P<d>[^"]+)"$', re.M),
        'date-released: "{d}"',
        lambda date: date,
    ),
    (
        "README.md",
        re.compile(r"^  year      = \{(?P<d>[^}]+)\},$", re.M),
        "  year      = {{{d}}},",
        lambda date: date[:4],  # BibTeX wants just the year
    ),
]


class ReleaseError(RuntimeError):
    """Something is wrong with the repo state; message is user-facing."""


def _read(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def _match_once(rel: str, pattern: re.Pattern, text: str) -> re.Match:
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ReleaseError(
            f"{rel}: expected exactly 1 match for {pattern.pattern!r}, found {len(matches)}. "
            "The file was reformatted — update the pattern in scripts/release.py."
        )
    return matches[0]


def current_release() -> tuple[str, str]:
    """Return (version, date) from the topmost CHANGELOG.md release heading."""
    history = released_versions()
    if not history:
        raise ReleaseError("CHANGELOG.md: no '## [version] - YYYY-MM-DD' heading found.")
    return history[0]


def released_versions() -> list[tuple[str, str]]:
    """Every (version, date) in CHANGELOG.md, newest first."""
    return [
        (m.group("v"), m.group("d")) for m in CHANGELOG_HEADING.finditer(_read("CHANGELOG.md"))
    ]


def previous_version(version: str) -> str | None:
    """The newest released version that isn't `version`.

    Not simply `history[0]`: the documented flow writes the CHANGELOG section
    *before* running `bump`, so by then the topmost heading already is the
    version being released. Comparing against it would yield an empty diff link.
    """
    for candidate, _ in released_versions():
        if candidate != version:
            return candidate
    return None


def audit(version: str, date: str) -> list[str]:
    """Return a list of human-readable drift reports (empty means in sync)."""
    problems = []
    for rel, pattern, _ in VERSION_SITES:
        found = _match_once(rel, pattern, _read(rel)).group("v")
        if found != version:
            problems.append(f"  {rel}: version is {found!r}, expected {version!r}")
    for rel, pattern, _, derive in DATE_SITES:
        found = _match_once(rel, pattern, _read(rel)).group("d")
        expected = derive(date)
        if found != expected:
            problems.append(f"  {rel}: date is {found!r}, expected {expected!r}")
    return problems


def unlinked_versions() -> list[str]:
    """Changelog versions with no `[x.y.z]: <url>` link reference at the bottom."""
    text = _read("CHANGELOG.md")
    declared = {m.group(1) for m in re.finditer(r"^\[([^\]]+)\]: http", text, re.M)}
    return [m.group("v") for m in CHANGELOG_HEADING.finditer(text) if m.group("v") not in declared]


def write_version(version: str, date: str, dry_run: bool) -> None:
    for rel, pattern, template in VERSION_SITES:
        text = _read(rel)
        match = _match_once(rel, pattern, text)
        updated = text[: match.start()] + template.format(v=version) + text[match.end() :]
        _emit(rel, match.group(0), template.format(v=version), updated, dry_run)
    for rel, pattern, template, derive in DATE_SITES:
        text = _read(rel)
        match = _match_once(rel, pattern, text)
        value = derive(date)
        updated = text[: match.start()] + template.format(d=value) + text[match.end() :]
        _emit(rel, match.group(0), template.format(d=value), updated, dry_run)


def _emit(rel: str, before: str, after: str, new_text: str, dry_run: bool) -> None:
    if before == after:
        print(f"  = {rel}: {after}")
        return
    print(f"  {'~' if dry_run else '+'} {rel}: {before}  ->  {after}")
    if not dry_run:
        (ROOT / rel).write_text(new_text, encoding="utf-8")


def ensure_changelog_section(version: str, date: str, dry_run: bool) -> None:
    """Insert a stub `## [version] - date` section if the author hasn't written one."""
    text = _read("CHANGELOG.md")
    if re.search(rf"^## \[{re.escape(version)}\] - ", text, re.M):
        return
    anchor = CHANGELOG_HEADING.search(text)
    insert_at = anchor.start() if anchor else len(text)
    stub = f"## [{version}] - {date}\n\n### Added\n\n- TODO: describe this release.\n\n"
    print(f"  {'~' if dry_run else '+'} CHANGELOG.md: inserted stub section for {version}")
    print("    !! Fill in the release notes before committing.")
    if not dry_run:
        (ROOT / "CHANGELOG.md").write_text(text[:insert_at] + stub + text[insert_at:], "utf-8")


def add_compare_link(version: str, previous: str, dry_run: bool) -> None:
    text = _read("CHANGELOG.md")
    if re.search(rf"^\[{re.escape(version)}\]: http", text, re.M):
        return
    link = f"[{version}]: {REPO_URL}/compare/v{previous}...v{version}\n"
    first_ref = re.search(r"^\[[^\]]+\]: http", text, re.M)
    if first_ref is None:
        text = text.rstrip("\n") + "\n\n" + link
    else:
        text = text[: first_ref.start()] + link + text[first_ref.start() :]
    print(f"  {'~' if dry_run else '+'} CHANGELOG.md: {link.strip()}")
    if not dry_run:
        (ROOT / "CHANGELOG.md").write_text(text, encoding="utf-8")


def unrelated_dirty_paths() -> list[str]:
    """Modified tracked files that a release is *not* expected to touch.

    Editing CHANGELOG.md before bumping is the documented flow, so an
    uncommitted changelog is fine; anything else means the bump would be mixed
    into unrelated work.
    """
    result = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:  # not a git checkout — nothing to guard
        return []
    managed = {"CHANGELOG.md", *(rel for rel, _, _ in VERSION_SITES)}
    dirty = (line[3:].strip() for line in result.stdout.splitlines() if line.strip())
    return sorted(path for path in dirty if path not in managed)


def cmd_check(_args: argparse.Namespace) -> int:
    version, date = current_release()
    print(f"CHANGELOG.md declares {version} ({date})")
    problems = audit(version, date)
    if problems:
        print("\nVersion strings are out of sync:")
        print("\n".join(problems))
        print("\nFix with: python scripts/release.py bump " + version)
        return 1
    print("All version strings are in sync.")
    missing = unlinked_versions()
    if missing:
        print(f"\nNote: no CHANGELOG compare link for: {', '.join(missing)}")
    return 0


def cmd_bump(args: argparse.Namespace) -> int:
    version = args.version.lstrip("v")
    if not SEMVER.match(version):
        raise ReleaseError(f"{version!r} is not a valid semantic version (expected X.Y.Z).")
    date = args.date or _dt.date.today().isoformat()
    if not ISO_DATE.match(date):
        raise ReleaseError(f"{date!r} is not an ISO date (expected YYYY-MM-DD).")

    previous = previous_version(version)
    if not args.dry_run:
        blocked = unrelated_dirty_paths()
        if blocked:
            raise ReleaseError(
                "These files have uncommitted changes unrelated to the release:\n"
                + "\n".join(f"  {path}" for path in blocked)
                + "\nCommit or stash them first, so the version bump lands as its own "
                "reviewable commit."
            )

    print(
        f"Bumping {previous or '(first release)'} -> {version} ({date})"
        + ("  [dry run]" if args.dry_run else "")
    )
    ensure_changelog_section(version, date, args.dry_run)
    if previous is not None:
        add_compare_link(version, previous, args.dry_run)
    write_version(version, date, args.dry_run)

    if args.dry_run:
        print("\nDry run: nothing written.")
        return 0

    print(
        f"\nDone. Remaining manual steps:\n"
        f"  1. Fill in / review the CHANGELOG.md section for {version}.\n"
        f"  2. pytest tests/ -v\n"
        f"  3. python scripts/release.py check\n"
        f"  4. git commit -am 'Release v{version}' && git tag v{version}\n"
        f"  5. git push && git push --tags\n"
        f"  6. python -m build && twine upload dist/*\n"
        f"  7. Confirm Zenodo minted a DOI for the new tag "
        f"(the concept DOI in CITATION.cff stays the same).\n"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("check", help="verify every version string matches CHANGELOG.md").set_defaults(
        func=cmd_check
    )

    bump = sub.add_parser("bump", help="set the version everywhere")
    bump.add_argument("version", help="new version, e.g. 0.1.16")
    bump.add_argument("--date", help="release date (YYYY-MM-DD); defaults to today")
    bump.add_argument("--dry-run", action="store_true", help="print changes without writing")
    bump.set_defaults(func=cmd_bump)

    args = parser.parse_args()
    try:
        return args.func(args)
    except ReleaseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
