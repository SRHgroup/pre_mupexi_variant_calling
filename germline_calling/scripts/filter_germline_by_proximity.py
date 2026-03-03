#!/usr/bin/env python3
"""Keep germline variants that are within N bp of any somatic (and optional RNA) site."""

from __future__ import annotations

import argparse
import gzip
from typing import Dict, Iterable, List, Set, TextIO, Tuple


def open_in(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def open_out(path: str) -> TextIO:
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w", encoding="utf-8")


def load_sites(path: str) -> Dict[str, List[int]]:
    sites: Dict[str, List[int]] = {}
    with open_in(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            sites.setdefault(chrom, []).append(pos)
    for chrom in sites:
        sites[chrom].sort()
    return sites


def has_nearby(sorted_positions: List[int], pos: int, distance: int) -> bool:
    # Linear scan is fine for moderate lists in per-sample usage.
    low = pos - distance
    high = pos + distance
    for p in sorted_positions:
        if p < low:
            continue
        if p > high:
            break
        return True
    return False


def filter_germline(germline_vcf: str, out_vcf: str, nearby_sites: Dict[str, List[int]], distance: int) -> Tuple[int, int]:
    total = 0
    kept = 0
    with open_in(germline_vcf) as fin, open_out(out_vcf) as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            total += 1
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            chrom_sites = nearby_sites.get(chrom)
            if chrom_sites is None:
                continue
            if has_nearby(chrom_sites, pos, distance):
                fout.write(line)
                kept += 1
    return total, kept


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter germline variants by adjacency to somatic/RNA variant positions")
    parser.add_argument("-g", "--germline", required=True, help="Input germline VCF(.gz)")
    parser.add_argument("-s", "--somatic", required=True, help="Input DNA somatic VCF(.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output filtered germline VCF(.gz)")
    parser.add_argument("-d", "--distance", type=int, default=33, help="Distance threshold in bp (default: 33)")
    parser.add_argument("--rna", default="", help="Optional RNA VCF(.gz) to include as adjacency source")
    args = parser.parse_args()

    nearby_sites = load_sites(args.somatic)
    if args.rna:
        rna_sites = load_sites(args.rna)
        for chrom, positions in rna_sites.items():
            nearby_sites.setdefault(chrom, []).extend(positions)
        for chrom in nearby_sites:
            nearby_sites[chrom].sort()

    total, kept = filter_germline(args.germline, args.output, nearby_sites, args.distance)
    print(f"[done] input_germline={total} kept={kept} distance={args.distance}")


if __name__ == "__main__":
    main()
