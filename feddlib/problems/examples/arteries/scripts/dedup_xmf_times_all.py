#!/usr/bin/env python3
"""Run XMF time deduplication across all subdirectories.

This is a convenience wrapper around tools/dedup_xmf_times.py.

Usage:
  python3 tools/dedup_xmf_times_all.py /path/to/root_dir

It will walk all subdirectories and, for every directory that contains at least
one of the target files (default: d_s.xmf, c.xmf), it will apply the same
"remove consecutive duplicate times" rewrite.

Backups:
  By default, creates .bak backups next to each rewritten XMF.

Examples:
  # Preview what would change
  python3 tools/dedup_xmf_times_all.py --dry-run /scratch/simulations

  # Actually fix everything
  python3 tools/dedup_xmf_times_all.py /scratch/simulations
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple


def _load_dedup_module(repo_root: Path):
    script = repo_root / "tools" / "dedup_xmf_times.py"
    if not script.exists():
        raise FileNotFoundError(f"Missing {script}")

    spec = importlib.util.spec_from_file_location("dedup_xmf_times", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {script}")
    mod = importlib.util.module_from_spec(spec)
    # Ensure the module is visible during execution (dataclasses inspects sys.modules).
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


def _find_candidate_dirs(root: Path, target_files: List[str]) -> List[Path]:
    """Return directories that contain at least one of target_files."""
    candidates: List[Path] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        filename_set = set(filenames)
        if any(name in filename_set for name in target_files):
            candidates.append(Path(dirpath))
    return candidates


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Deduplicate duplicate XMF times in all subdirectories (d_s.xmf, c.xmf by default).")
    parser.add_argument("root_dir", type=Path, help="Root directory to search recursively")
    parser.add_argument("--files", nargs="*", default=["d_s.xmf", "c.xmf"], help="XMF filenames to process in each directory")
    parser.add_argument("--no-backup", action="store_true", help="Do not create .bak backups")
    parser.add_argument("--dry-run", action="store_true", help="Analyze and report changes but do not modify files")

    args = parser.parse_args(argv)

    root_dir: Path = args.root_dir
    if not root_dir.is_dir():
        print(f"error: not a directory: {root_dir}", file=sys.stderr)
        return 2

    # Assume this script lives in <repo>/tools; repo root is parent of tools.
    repo_root = Path(__file__).resolve().parent.parent
    dedup = _load_dedup_module(repo_root)

    candidate_dirs = _find_candidate_dirs(root_dir, list(args.files))
    if not candidate_dirs:
        print(f"No directories under {root_dir} contained any of: {', '.join(args.files)}")
        return 0

    total_dirs = 0
    total_files = 0
    total_removed = 0
    rc = 0

    for d in sorted(candidate_dirs):
        total_dirs += 1
        changed_any = False
        for name in args.files:
            path = d / name
            if not path.exists():
                continue
            try:
                num_blocks, removed = dedup.process_file(path, make_backup=(not args.no_backup), dry_run=args.dry_run)
                total_files += 1
                total_removed += removed
                if removed:
                    changed_any = True
                    action = "would fix" if args.dry_run else "fixed"
                    print(f"{action}: {path} (timesteps={num_blocks}, removed_duplicates={removed})")
            except Exception as e:
                rc = 1
                print(f"error: {path}: {e}", file=sys.stderr)

        if changed_any:
            pass

    summary_action = "Would update" if args.dry_run else "Updated"
    print(f"{summary_action} {total_files} file(s) across {total_dirs} directorie(s); removed {total_removed} duplicate timestep(s).")

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
