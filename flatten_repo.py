#!/usr/bin/env python3
"""
flatten_repo.py
---------------
Copy every file in a repository into ONE flat folder (no subfolders), and
write a Markdown file that records the original directory structure plus a
map from each flattened filename back to its original path.

Run from the repo root:

    python flatten_repo.py                 # -> ./_flat_export/
    python flatten_repo.py --out upload     # custom output folder
    python flatten_repo.py --clean          # wipe output folder first
    python flatten_repo.py --root ../myrepo # flatten a different repo

Filename collisions (e.g. many __init__.py) are handled automatically:
a file keeps its plain name if that name is unique, otherwise it is
renamed using its full relative path (src/models/foo.py -> src__models__foo.py),
so nothing is ever silently overwritten. Collisions are checked
case-insensitively, which matters on Windows / OneDrive filesystems.
"""

from __future__ import annotations

import argparse
import os
import shutil
from collections import Counter
from pathlib import Path

OUTPUT_DIRNAME_DEFAULT = "_flat_export"
STRUCTURE_FILENAME = "_STRUCTURE.md"

# Directories that are never worth copying.
SKIP_DIRS = {
    ".git", ".hg", ".svn",
    "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".tox",
    ".venv", "venv", "env", ".env",
    "node_modules", ".idea", ".vscode", ".ipynb_checkpoints",
    "build", "dist", ".eggs", ".cache","tests","examples",
    ".jax_cache", ".pytest_cache", ".ruff_cache", ".docs", "tests", "Smoothed_TEfixed_linear", "build"
}

# Individual files to skip.
SKIP_SUFFIXES = {".pyc", ".pyo", ".pyd", ".so", ".o", ".a", ".class"}
SKIP_NAMES = {".DS_Store", "Thumbs.db"}


def should_skip_file(path: Path) -> bool:
    if path.name in SKIP_NAMES:
        return True
    if path.suffix.lower() in SKIP_SUFFIXES:
        return True
    return False


def build_tree(directory: Path, output_dirname: str, prefix: str = "") -> list[str]:
    """Return an ASCII tree of `directory`, skipping junk and the output folder."""
    try:
        children = list(directory.iterdir())
    except PermissionError:
        return []

    entries = [
        p for p in children
        if not (p.is_dir() and (p.name in SKIP_DIRS or p.name == output_dirname))
        and not (p.is_file() and should_skip_file(p))
    ]
    # Directories first, then files, each alphabetically (case-insensitive).
    entries.sort(key=lambda p: (p.is_file(), p.name.lower()))

    lines: list[str] = []
    for i, p in enumerate(entries):
        last = i == len(entries) - 1
        connector = "└── " if last else "├── "
        lines.append(prefix + connector + p.name + ("/" if p.is_dir() else ""))
        if p.is_dir():
            extension = "    " if last else "│   "
            lines.extend(build_tree(p, output_dirname, prefix + extension))
    return lines


def collect_files(root: Path, output_dirname: str) -> list[Path]:
    """Walk `root`, pruning skip-dirs and the output folder, return file paths."""
    files: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        # Prune in place so os.walk never descends into these.
        dirnames[:] = [
            d for d in dirnames
            if d not in SKIP_DIRS and d != output_dirname
        ]
        for fn in filenames:
            p = Path(dirpath) / fn
            if not should_skip_file(p):
                files.append(p)
    return files


def assign_flat_names(root: Path, files: list[Path]) -> dict[Path, str]:
    """Map each source file to a unique flat filename."""
    rels = {f: f.relative_to(root) for f in files}
    basename_counts = Counter(r.name.lower() for r in rels.values())

    mapping: dict[Path, str] = {}
    used: set[str] = set()

    # Stable order: by relative path.
    for f in sorted(files, key=lambda p: str(rels[p]).lower()):
        rel = rels[f]
        if basename_counts[rel.name.lower()] == 1:
            name = rel.name
        else:
            # Encode the full relative path into the filename.
            name = "__".join(rel.parts)

        # Final safety net against any residual clash.
        base, n = name, 1
        while name.lower() in used:
            stem, suffix = Path(base).stem, Path(base).suffix
            name = f"{stem}_{n}{suffix}"
            n += 1
        used.add(name.lower())
        mapping[f] = name
    return mapping


def write_structure_file(
    dest: Path, root: Path, output_dirname: str, mapping: dict[Path, str]
) -> None:
    tree_lines = build_tree(root, output_dirname)
    repo_name = root.resolve().name or "repo"

    # Column-aligned flat-name -> original-path map.
    rows = sorted(
        ((name, str(f.relative_to(root)).replace(os.sep, "/")) for f, name in mapping.items()),
        key=lambda r: r[1].lower(),
    )
    width = max((len(name) for name, _ in rows), default=0)
    map_lines = [f"{name.ljust(width)}  <-  {orig}" for name, orig in rows]

    content = (
        f"# Original structure of `{repo_name}`\n\n"
        f"{len(mapping)} files flattened into a single folder.\n\n"
        "## Directory tree\n\n"
        "```\n"
        f"{repo_name}/\n"
        + "\n".join(tree_lines)
        + "\n```\n\n"
        "## Flattened file map\n\n"
        "Each flattened filename and where it came from:\n\n"
        "```\n"
        + "\n".join(map_lines)
        + "\n```\n"
    )
    dest.write_text(content, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default=".", help="repo root to flatten (default: .)")
    parser.add_argument("--out", default=OUTPUT_DIRNAME_DEFAULT,
                        help=f"output folder name (default: {OUTPUT_DIRNAME_DEFAULT})")
    parser.add_argument("--clean", action="store_true",
                        help="delete the output folder before copying")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    output_dirname = args.out
    out_dir = root / output_dirname

    if not root.is_dir():
        raise SystemExit(f"error: {root} is not a directory")

    if args.clean and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = collect_files(root, output_dirname)
    if not files:
        raise SystemExit("no files found to copy (everything was skipped?)")

    mapping = assign_flat_names(root, files)

    renamed = 0
    for src, name in mapping.items():
        shutil.copy2(src, out_dir / name)
        if name != src.name:
            renamed += 1

    write_structure_file(out_dir / STRUCTURE_FILENAME, root, output_dirname, mapping)

    print(f"Copied {len(files)} files into {out_dir}")
    print(f"  {renamed} renamed to avoid collisions")
    print(f"  structure written to {out_dir / STRUCTURE_FILENAME}")


if __name__ == "__main__":
    main()