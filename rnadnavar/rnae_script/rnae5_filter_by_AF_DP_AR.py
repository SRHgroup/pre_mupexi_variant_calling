#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, sys
from typing import Dict, List, Optional, TextIO, Tuple

def open_in(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def open_out(path: str) -> TextIO:
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

def parse_int(x: str) -> Optional[int]:
    if x is None or x == "." or x == "":
        return None
    try:
        return int(x)
    except ValueError:
        return None

def parse_ad(ad_str: str) -> Optional[List[int]]:
    if ad_str is None or ad_str == "." or ad_str == "":
        return None
    parts = ad_str.split(",")
    out: List[int] = []
    for p in parts:
        v = parse_int(p)
        if v is None:
            return None
        out.append(v)
    return out

def parse_sample(fmt_keys: List[str], sample_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    fields = sample_str.split(":")
    for k, v in zip(fmt_keys, fields):
        d[k] = v
    return d

def dp_alt_vaf(sample: Dict[str, str]) -> Tuple[Optional[int], Optional[int], Optional[float]]:
    """
    Compute DP, ALT_READS, VAF from FORMAT fields (prefer AD; DP fallback to sum(AD)).
    VAF uses AD: ALT / (REF+ALT). If AD missing -> returns (DP, None, None).
    """
    dp = parse_int(sample.get("DP", "."))
    ad = parse_ad(sample.get("AD", "."))
    if ad is None or len(ad) < 2:
        return dp, None, None

    ref = ad[0]
    alt = sum(ad[1:])  # supports multi-ALT if present
    denom = ref + alt
    vaf = (alt / denom) if denom > 0 else None

    if dp is None:
        dp = denom
    return dp, alt, vaf

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output VCF (.vcf or .vcf.gz)")

    ap.add_argument("--normal-sample", default="DNA_NORMAL")
    ap.add_argument("--rna-sample", default="RNA_TUMOUR")

    ap.add_argument("--normal-min-dp", type=int, default=10)
    ap.add_argument("--normal-max-vaf", type=float, default=0.05)  # strict <
    ap.add_argument("--normal-max-alt", type=int, default=1)

    ap.add_argument("--rna-min-dp", type=int, default=10)
    ap.add_argument("--rna-min-vaf", type=float, default=0.05)     # >=
    ap.add_argument("--rna-min-alt", type=int, default=3)

    ap.add_argument("--stats", action="store_true", help="Print filtering stats to stderr")
    args = ap.parse_args()

    normal_idx: Optional[int] = None
    rna_idx: Optional[int] = None

    n_total = 0
    n_kept = 0
    drop_missing = 0
    drop_normal = 0
    drop_rna = 0

    with open_in(args.input) as fin, open_out(args.output) as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                for i, s in enumerate(samples):
                    if s == args.normal_sample:
                        normal_idx = i
                    if s == args.rna_sample:
                        rna_idx = i
                if normal_idx is None or rna_idx is None:
                    raise SystemExit(
                        f"ERROR: Could not find required samples. "
                        f"normal='{args.normal_sample}' rna='{args.rna_sample}'"
                    )
                fout.write(line)
                continue

            if not line or line[0] == "#":
                fout.write(line)
                continue

            toks = line.rstrip("\n").split("\t")
            if len(toks) < 10:
                continue

            n_total += 1

            fmt_keys = toks[8].split(":")
            sample_fields = toks[9:]

            normal = parse_sample(fmt_keys, sample_fields[normal_idx])
            rna = parse_sample(fmt_keys, sample_fields[rna_idx])

            n_dp, n_alt, n_vaf = dp_alt_vaf(normal)
            r_dp, r_alt, r_vaf = dp_alt_vaf(rna)

            # Require AD-derived alt/vaf in both, and dp present/computable
            if n_dp is None or r_dp is None or n_alt is None or r_alt is None or n_vaf is None or r_vaf is None:
                drop_missing += 1
                continue

            # Normal criteria
            if not (n_dp >= args.normal_min_dp and n_alt <= args.normal_max_alt and n_vaf < args.normal_max_vaf):
                drop_normal += 1
                continue

            # RNA criteria
            if not (r_dp >= args.rna_min_dp and r_alt >= args.rna_min_alt and r_vaf >= args.rna_min_vaf):
                drop_rna += 1
                continue

            fout.write(line)
            n_kept += 1

    if args.stats:
        print("filter_stats", file=sys.stderr)
        print(f"  input_total\t{n_total}", file=sys.stderr)
        print(f"  kept\t{n_kept}", file=sys.stderr)
        print(f"  dropped_missing_AD_or_DP\t{drop_missing}", file=sys.stderr)
        print(f"  dropped_normal_fail\t{drop_normal}", file=sys.stderr)
        print(f"  dropped_rna_fail\t{drop_rna}", file=sys.stderr)

if __name__ == "__main__":
    main()
