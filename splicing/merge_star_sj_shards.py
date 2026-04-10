#!/usr/bin/env python3
"""Merge STAR SJ.out.tab shard files into one STAR-style file per sample."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


SHARD_RE = re.compile(r"^(?P<base>.+)\.(?P<shard>\d{4})\.SJ\.out\.tab(?P<gz>\.gz)?$")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Recursively find STAR SJ.out.tab shard files such as "
            "Pat21_RNA_TUMOUR.0001.SJ.out.tab and merge them into a single "
            "STAR-style output named Pat21_RNA_TUMOUR.SJ.out.tab."
        )
    )
    ap.add_argument("--root", required=True, help="Root directory to scan recursively")
    ap.add_argument(
        "--sample-filter",
        action="append",
        default=[],
        help=(
            "Substring filter applied to the output path/base name. "
            "Repeatable; any match is kept."
        ),
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing merged outputs",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Report planned merges without writing files",
    )
    return ap.parse_args()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def chrom_sort_key(chrom: str) -> Tuple[int, object]:
    c = chrom
    if c.startswith("chr"):
        c = c[3:]
    if c.isdigit():
        return (0, int(c))
    if c == "X":
        return (0, 23)
    if c == "Y":
        return (0, 24)
    if c in {"M", "MT"}:
        return (0, 25)
    return (1, chrom)


def sample_matches(group_out: Path, filters: Sequence[str]) -> bool:
    if not filters:
        return True
    haystacks = {group_out.name.lower(), str(group_out).lower(), group_out.stem.lower()}
    return any(f.lower() in text for f in filters for text in haystacks)


def discover_groups(root: Path, filters: Sequence[str]) -> Dict[Tuple[Path, str], List[Tuple[int, Path]]]:
    groups: Dict[Tuple[Path, str], List[Tuple[int, Path]]] = defaultdict(list)
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        match = SHARD_RE.match(path.name)
        if not match:
            continue
        base = match.group("base")
        shard_no = int(match.group("shard"))
        groups[(path.parent, base)].append((shard_no, path))

    filtered: Dict[Tuple[Path, str], List[Tuple[int, Path]]] = {}
    for key, shard_paths in groups.items():
        parent, base = key
        all_gz = all(str(p).endswith(".gz") for _, p in shard_paths)
        out_name = f"{base}.SJ.out.tab" + (".gz" if all_gz else "")
        out_path = parent / out_name
        if sample_matches(out_path, filters):
            filtered[key] = sorted(shard_paths, key=lambda item: item[0])
    return filtered


def merge_group(shard_paths: Sequence[Path]) -> List[List[str]]:
    merged: Dict[Tuple[str, str, str, str, str, str], List[object]] = {}
    for shard_path in shard_paths:
        with open_text(shard_path, "rt") as fh:
            for lineno, raw_line in enumerate(fh, start=1):
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                cols = line.split("\t")
                if len(cols) != 9:
                    raise SystemExit(
                        f"ERROR: expected 9 tab-delimited columns in {shard_path}:{lineno}, got {len(cols)}"
                    )
                key = tuple(cols[:6])
                try:
                    uniq = int(cols[6])
                    multi = int(cols[7])
                    overhang = int(cols[8])
                except ValueError as exc:
                    raise SystemExit(
                        f"ERROR: non-integer count field in {shard_path}:{lineno}: {exc}"
                    ) from exc

                rec = merged.get(key)
                if rec is None:
                    merged[key] = [list(cols[:6]), uniq, multi, overhang]
                else:
                    rec[1] += uniq
                    rec[2] += multi
                    rec[3] = max(rec[3], overhang)

    ordered = sorted(
        merged.values(),
        key=lambda rec: (
            chrom_sort_key(rec[0][0]),
            int(rec[0][1]),
            int(rec[0][2]),
            rec[0][3],
            rec[0][4],
            rec[0][5],
        ),
    )
    return [
        rec[0] + [str(rec[1]), str(rec[2]), str(rec[3])]
        for rec in ordered
    ]


def write_rows(path: Path, rows: Iterable[Sequence[str]]) -> None:
    with open_text(path, "wt") as fh:
        for row in rows:
            fh.write("\t".join(row))
            fh.write("\n")


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"ERROR: root does not exist: {root}")
    if not root.is_dir():
        raise SystemExit(f"ERROR: root is not a directory: {root}")

    groups = discover_groups(root, args.sample_filter)
    if not groups:
        print(f"[info] no STAR SJ shard groups found under {root}", file=sys.stderr)
        return 0

    merged_groups = 0
    skipped_groups = 0
    for (parent, base), shard_items in sorted(groups.items(), key=lambda item: (str(item[0][0]), item[0][1])):
        shard_paths = [path for _, path in shard_items]
        all_gz = all(str(path).endswith(".gz") for path in shard_paths)
        out_path = parent / (f"{base}.SJ.out.tab" + (".gz" if all_gz else ""))

        if out_path.exists() and not args.force:
            print(f"[skip] output exists: {out_path}", file=sys.stderr)
            skipped_groups += 1
            continue

        print(
            f"[merge] {out_path} <- {len(shard_paths)} shard(s)",
            file=sys.stderr,
        )
        for shard_path in shard_paths:
            print(f"  [in] {shard_path}", file=sys.stderr)

        if args.dry_run:
            merged_groups += 1
            continue

        rows = merge_group(shard_paths)
        write_rows(out_path, rows)
        print(f"  [out] {out_path} ({len(rows)} junctions)", file=sys.stderr)
        merged_groups += 1

    print(
        f"[done] merged_groups={merged_groups} skipped_groups={skipped_groups}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
