#!/usr/bin/env python3
"""Deduplicate consecutive duplicate time steps in FEDDLib ParaView .xmf outputs.

Why this exists:
Some runs can end up writing two consecutive <Grid> blocks with the same printed
Time value (e.g. "1500" twice). ParaView may mis-handle duplicate timestamps.

This tool rewrites the XMF by keeping only the *last* occurrence of any
consecutive duplicate time value.

Default behavior targets the common FEDDLib outputs: d_s.xmf and c.xmf.

Usage:
  python3 tools/dedup_xmf_times.py /path/to/results_dir

This will:
  - create backups: d_s.xmf.bak, c.xmf.bak
  - write fixed files in-place

Exit status is non-zero if any requested file exists but cannot be processed.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


_TIME_RE = re.compile(r"TimeType=\"Single\"\s+Value=\"([^\"]+)\"")
_GRID_OPEN_RE = re.compile(r"<Grid\b")
_GRID_CLOSE_RE = re.compile(r"</Grid>")


@dataclass(frozen=True)
class GridBlock:
    time_str: str
    time_val: float
    text: str


def _parse_float(s: str) -> float:
    try:
        return float(s)
    except Exception:
        # If it can't be parsed, treat as NaN-ish; still dedup by string.
        return float("nan")


def _is_close(a: float, b: float, *, rel: float = 1e-12, abs_: float = 0.0) -> bool:
    # Don’t import math.isclose to avoid any surprises with NaN.
    if a != a or b != b:  # NaN
        return False
    diff = abs(a - b)
    if diff <= abs_:
        return True
    return diff <= rel * max(abs(a), abs(b), 1.0)


def _extract_inner_grids(xmf_text: str) -> Tuple[str, List[GridBlock], str]:
    """Split XMF into header, list of inner Grid blocks, and footer.

    Assumptions:
      - The file contains one outer temporal collection Grid.
      - Inner per-timestep grids are direct children of that collection.

    We locate inner grids by finding "<Grid Name=\"Mesh..." (most FEDDLib files
    use that pattern) and then capturing until its matching </Grid>.

    If pattern search fails, we fall back to a depth-based scan capturing any
    Grid blocks at depth==2.
    """

    # Fast path: explicit “Mesh” named grids.
    mesh_start = re.search(r"<Grid\s+Name=\"Mesh", xmf_text)
    if mesh_start:
        header = xmf_text[: mesh_start.start()]
        remainder = xmf_text[mesh_start.start() :]

        blocks: List[GridBlock] = []
        idx = 0
        while True:
            m = re.search(r"<Grid\s+Name=\"Mesh", remainder[idx:])
            if not m:
                break
            start = idx + m.start()

            # Capture a complete <Grid> ... </Grid> block by counting Grid tags.
            local_depth = 0
            end = None
            for mm in re.finditer(r"<Grid\b|</Grid>", remainder[start:]):
                token = mm.group(0)
                if token.startswith("<Grid"):
                    local_depth += 1
                else:
                    local_depth -= 1
                    if local_depth == 0:
                        end = start + mm.end()
                        break
            if end is None:
                raise ValueError("Unterminated <Grid> block while parsing XMF")

            block_text = remainder[start:end]
            time_match = _TIME_RE.search(block_text)
            if not time_match:
                raise ValueError("Could not find <Time ... Value=...> inside a Mesh grid")
            time_str = time_match.group(1)
            blocks.append(GridBlock(time_str=time_str, time_val=_parse_float(time_str), text=block_text))

            idx = end

        footer = remainder[idx:]
        return header, blocks, footer

    # Fallback: depth-based scan.
    lines = xmf_text.splitlines(keepends=True)
    header_parts: List[str] = []
    footer_parts: List[str] = []
    blocks: List[GridBlock] = []

    outer_depth = 0
    capturing = False
    cur_block_lines: List[str] = []
    cur_local_depth = 0

    # We consider inner blocks those that start when overall depth becomes 2.
    # The first Grid (outer collection) brings depth to 1.
    for line in lines:
        opens = len(_GRID_OPEN_RE.findall(line))
        closes = len(_GRID_CLOSE_RE.findall(line))

        if not capturing:
            # Detect start of an inner grid when already inside outer collection.
            if outer_depth == 1 and "<Grid" in line and opens > closes:
                capturing = True
                cur_block_lines = [line]
                cur_local_depth = opens - closes
            else:
                header_parts.append(line)
        else:
            cur_block_lines.append(line)
            cur_local_depth += opens - closes
            if cur_local_depth == 0:
                capturing = False
                block_text = "".join(cur_block_lines)
                time_match = _TIME_RE.search(block_text)
                if not time_match:
                    raise ValueError("Could not find <Time ... Value=...> inside an inner grid")
                time_str = time_match.group(1)
                blocks.append(GridBlock(time_str=time_str, time_val=_parse_float(time_str), text=block_text))

        outer_depth += opens - closes

    if capturing:
        raise ValueError("Unterminated <Grid> block while parsing XMF")

    # In fallback mode, everything went to header_parts; we don’t have a clean footer.
    # Treat footer as empty.
    return "".join(header_parts), blocks, "".join(footer_parts)


def dedup_consecutive_times(blocks: List[GridBlock]) -> Tuple[List[GridBlock], int]:
    """Keep the last of any consecutive duplicate times."""
    if not blocks:
        return blocks, 0

    out: List[GridBlock] = []
    removed = 0

    for block in blocks:
        if not out:
            out.append(block)
            continue

        prev = out[-1]
        # Primary: exact string equality (most relevant for ParaView duplicate time display)
        # Secondary: numeric closeness to catch e.g. "1500" vs "1500.0".
        if block.time_str == prev.time_str or _is_close(block.time_val, prev.time_val, rel=0.0, abs_=0.0):
            out[-1] = block
            removed += 1
        else:
            out.append(block)

    return out, removed


def process_file(path: Path, *, make_backup: bool = True, dry_run: bool = False) -> Tuple[int, int]:
    """Return (num_blocks, num_removed)."""
    text = path.read_text(encoding="utf-8", errors="replace")
    header, blocks, footer = _extract_inner_grids(text)

    deduped, removed = dedup_consecutive_times(blocks)

    if removed == 0:
        return len(blocks), 0

    new_text = header + "".join(b.text for b in deduped) + footer

    if dry_run:
        return len(blocks), removed

    if make_backup:
        backup = path.with_suffix(path.suffix + ".bak")
        if not backup.exists():
            shutil.copy2(path, backup)

    # Atomic-ish write.
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as tmp:
        tmp.write(new_text)
        tmp_path = Path(tmp.name)

    tmp_path.replace(path)

    return len(blocks), removed


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Remove consecutive duplicate time entries from FEDDLib XMF files (d_s.xmf, c.xmf).")
    parser.add_argument("results_dir", type=Path, help="Directory containing d_s.xmf and c.xmf")
    parser.add_argument("--files", nargs="*", default=["d_s.xmf", "c.xmf"], help="XMF filenames to process (default: d_s.xmf c.xmf)")
    parser.add_argument("--no-backup", action="store_true", help="Do not create .bak backups")
    parser.add_argument("--dry-run", action="store_true", help="Analyze and report changes but do not modify files")

    args = parser.parse_args(argv)

    results_dir: Path = args.results_dir
    if not results_dir.is_dir():
        print(f"error: not a directory: {results_dir}", file=sys.stderr)
        return 2

    rc = 0
    for name in args.files:
        path = results_dir / name
        if not path.exists():
            print(f"skip: {path} (not found)")
            continue
        try:
            num_blocks, removed = process_file(path, make_backup=(not args.no_backup), dry_run=args.dry_run)
            if removed:
                action = "would fix" if args.dry_run else "fixed"
                print(f"{action}: {path} (timesteps={num_blocks}, removed_duplicates={removed})")
            else:
                print(f"ok: {path} (timesteps={num_blocks}, no duplicates)")
        except Exception as e:
            print(f"error: {path}: {e}", file=sys.stderr)
            rc = 1

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
