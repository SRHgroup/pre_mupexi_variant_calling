#!/usr/bin/env python3
"""rnae4_summarise_RNA_metrics.py

Collapse multiple RNA tumour/tumor samples (e.g. *_RNA_TUMOR.0001..0004) into ONE RNA sample,
and keep exactly two samples total in the output: DNA normal and RNA tumour/tumor.

This version is *label-aware*:
- It can detect prefixed sample names by suffix match (e.g. 01-CH-L_01-CH-L_DNA_NORMAL).
- It supports both spellings (TUMOR/TUMOUR) by default.

Aggregation (per-variant, across all RNA tumour samples):
- AD, DP, F1R2, F2R1, FAD, SB: summed element-wise where present
- AF: recomputed from summed AD as ALT_SUM / DP_SUM
- GT: 0/1 if ALT_SUM>0 else 0/0 (./. if DP_SUM==0 or missing)

All other FORMAT keys are dropped from BOTH samples.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from typing import Dict, List, Optional, TextIO

WANTED_KEYS_ORDER = ["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"]


def open_in(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def open_out(path: str) -> TextIO:
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def safe_int(x: str) -> Optional[int]:
    if x is None or x == "." or x == "":
        return None
    try:
        return int(x)
    except ValueError:
        return None


def parse_int_list(val: str) -> Optional[List[int]]:
    if val is None or val == "." or val == "":
        return None
    parts = val.split(",")
    out: List[int] = []
    for p in parts:
        i = safe_int(p)
        if i is None:
            return None
        out.append(i)
    return out


def sum_lists(a: Optional[List[int]], b: Optional[List[int]]) -> Optional[List[int]]:
    if a is None:
        return b
    if b is None:
        return a
    if len(a) != len(b):
        return None
    return [x + y for x, y in zip(a, b)]


def parse_sample(fmt_keys: List[str], sample_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    if sample_str in (".", "./.", ".:."):
        return d
    fields = sample_str.split(":")
    for k, v in zip(fmt_keys, fields):
        d[k] = v
    return d


def recompute_af_from_ad_dp(ad: Optional[List[int]], dp: Optional[int]) -> Optional[float]:
    if ad is None or len(ad) < 2:
        return None
    if dp is None:
        dp = sum(ad)
    if dp == 0:
        return None
    alt_sum = sum(ad[1:])
    return alt_sum / dp


def fmt_float(x: Optional[float]) -> str:
    if x is None:
        return "."
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    return s if s else "0"


def build_sample_string(keys: List[str], values: Dict[str, str]) -> str:
    return ":".join(values.get(k, ".") for k in keys)


def tumor_spellings(label: str) -> List[str]:
    """Return label plus the alternate spelling if it contains TUMOR/TUMOUR."""
    labels = [label]
    if "TUMOR" in label and "TUMOUR" not in labels:
        labels.append(label.replace("TUMOR", "TUMOUR"))
    if "TUMOUR" in label and "TUMOR" not in labels:
        labels.append(label.replace("TUMOUR", "TUMOR"))
    # de-dup preserving order
    out: List[str] = []
    for x in labels:
        if x not in out:
            out.append(x)
    return out


def find_sample_index_by_suffix(samples: List[str], suffixes: List[str]) -> Optional[int]:
    """Find first sample that equals a suffix or ends with _suffix."""
    for i, s in enumerate(samples):
        for suf in suffixes:
            if s == suf or s.endswith("_" + suf):
                return i
    return None


def find_rna_indices(samples: List[str], rna_suffixes: List[str]) -> List[int]:
    """Match RNA samples by suffix, allowing optional .0001 lane suffix."""
    idxs: List[int] = []
    # match: <anything>_<RNA_TUMOR> or <RNA_TUMOR> optionally followed by .dddd
    pats = [re.compile(rf"(?:^|_){re.escape(suf)}(?:\\.\\d+)?$") for suf in rna_suffixes]
    for i, s in enumerate(samples):
        if any(p.search(s) for p in pats):
            idxs.append(i)
    return idxs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input multi-sample VCF (.vcf or .vcf.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output 2-sample VCF (.vcf or .vcf.gz)")

    # Detection: either explicit sample names OR label-based detection
    ap.add_argument("--normal-sample", default=None, help="Exact DNA normal sample name in input header")
    ap.add_argument("--rna-sample", default=None, help="Exact RNA sample name/prefix in input header (optional)")

    ap.add_argument("--normal-label", default="DNA_NORMAL", help="Suffix label to detect DNA normal sample")
    ap.add_argument("--rna-label", default="RNA_TUMOR", help="Suffix label to detect RNA tumour/tumor samples")

    # Standardize output header
    ap.add_argument("--out-normal-name", default=None, help="Output name for DNA normal sample")
    ap.add_argument("--out-rna-name", default=None, help="Output name for collapsed RNA sample")

    args = ap.parse_args()

    normal_label_opts = [args.normal_label]
    rna_label_opts = tumor_spellings(args.rna_label)

    if args.out_normal_name is None:
        args.out_normal_name = args.normal_label
    if args.out_rna_name is None:
        args.out_rna_name = args.rna_label

    wrote_header = False
    normal_idx: Optional[int] = None
    rna_indices: List[int] = []

    n_in = 0
    n_out = 0

    with open_in(args.input) as fin, open_out(args.output) as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                fixed = cols[:9]
                samples = cols[9:]

                # locate DNA normal
                if args.normal_sample:
                    try:
                        normal_idx = samples.index(args.normal_sample)
                    except ValueError:
                        raise SystemExit(
                            f"ERROR: Could not find normal sample named '{args.normal_sample}' in header."
                        )
                else:
                    normal_idx = find_sample_index_by_suffix(samples, normal_label_opts)
                    if normal_idx is None:
                        raise SystemExit(
                            f"ERROR: Could not detect normal sample by label(s) {normal_label_opts} in header."
                        )

                # locate RNA samples
                if args.rna_sample:
                    # if exact sample provided, collapse just that sample
                    if args.rna_sample in samples:
                        rna_indices = [samples.index(args.rna_sample)]
                    else:
                        # treat as prefix (strip optional .0001): collapse all that start with it
                        prefix = args.rna_sample.split(".")[0]
                        rna_indices = [i for i, s in enumerate(samples) if s.split(".")[0] == prefix]
                else:
                    rna_indices = find_rna_indices(samples, rna_label_opts)

                # remove normal from RNA indices if it somehow matched
                rna_indices = [i for i in rna_indices if i != normal_idx]

                if not rna_indices:
                    raise SystemExit(
                        f"ERROR: No RNA samples detected using rna-label(s) {rna_label_opts} (or --rna-sample)."
                    )

                fout.write("\t".join(fixed + [args.out_normal_name, args.out_rna_name]) + "\n")
                wrote_header = True
                continue

            if not wrote_header:
                fout.write(line)
                continue

            if not line.strip() or line[0] == "#":
                continue

            n_in += 1
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 10:
                continue

            fmt_keys = toks[8].split(":")
            sample_fields = toks[9:]

            dn = parse_sample(fmt_keys, sample_fields[normal_idx])

            # Aggregate RNA
            agg_ad: Optional[List[int]] = None
            agg_f1r2: Optional[List[int]] = None
            agg_f2r1: Optional[List[int]] = None
            agg_fad: Optional[List[int]] = None
            agg_sb: Optional[List[int]] = None
            dp_sum: Optional[int] = 0
            dp_any = False

            for ridx in rna_indices:
                rs = parse_sample(fmt_keys, sample_fields[ridx])

                agg_ad = sum_lists(agg_ad, parse_int_list(rs.get("AD", ".")))
                agg_f1r2 = sum_lists(agg_f1r2, parse_int_list(rs.get("F1R2", ".")))
                agg_f2r1 = sum_lists(agg_f2r1, parse_int_list(rs.get("F2R1", ".")))
                agg_fad = sum_lists(agg_fad, parse_int_list(rs.get("FAD", ".")))
                agg_sb = sum_lists(agg_sb, parse_int_list(rs.get("SB", ".")))

                dp = safe_int(rs.get("DP", "."))
                if dp is not None:
                    dp_sum = (dp_sum or 0) + dp
                    dp_any = True

            dp_sum = dp_sum if dp_any else None

            if dp_sum is None and agg_ad is not None:
                dp_sum = sum(agg_ad)

            rna_af = recompute_af_from_ad_dp(agg_ad, dp_sum)

            if dp_sum is None or dp_sum == 0 or agg_ad is None:
                rna_gt = "./."
            else:
                alt_sum = sum(agg_ad[1:]) if len(agg_ad) > 1 else 0
                rna_gt = "0/1" if alt_sum > 0 else "0/0"

            keep_keys = [k for k in WANTED_KEYS_ORDER if k in fmt_keys]

            # DNA normal output
            dn_out: Dict[str, str] = {}
            dn_ad = parse_int_list(dn.get("AD", ".")) if "AD" in keep_keys else None
            dn_dp = safe_int(dn.get("DP", ".")) if "DP" in keep_keys else None
            if dn_dp is None and dn_ad is not None:
                dn_dp = sum(dn_ad)
            dn_af = recompute_af_from_ad_dp(dn_ad, dn_dp) if "AF" in keep_keys else None

            for k in keep_keys:
                if k == "GT":
                    dn_out[k] = dn.get("GT", ".")
                elif k == "AD":
                    dn_out[k] = dn.get("AD", ".")
                elif k == "DP":
                    dn_out[k] = str(dn_dp) if dn_dp is not None else dn.get("DP", ".")
                elif k == "AF":
                    dn_out[k] = fmt_float(dn_af)
                else:
                    dn_out[k] = dn.get(k, ".")

            # RNA collapsed output
            rn_out: Dict[str, str] = {}
            for k in keep_keys:
                if k == "GT":
                    rn_out[k] = rna_gt
                elif k == "AD":
                    rn_out[k] = ",".join(map(str, agg_ad)) if agg_ad is not None else "."
                elif k == "DP":
                    rn_out[k] = str(dp_sum) if dp_sum is not None else "."
                elif k == "AF":
                    rn_out[k] = fmt_float(rna_af)
                elif k == "F1R2":
                    rn_out[k] = ",".join(map(str, agg_f1r2)) if agg_f1r2 is not None else "."
                elif k == "F2R1":
                    rn_out[k] = ",".join(map(str, agg_f2r1)) if agg_f2r1 is not None else "."
                elif k == "FAD":
                    rn_out[k] = ",".join(map(str, agg_fad)) if agg_fad is not None else "."
                elif k == "SB":
                    rn_out[k] = ",".join(map(str, agg_sb)) if agg_sb is not None else "."
                else:
                    rn_out[k] = "."

            toks[8] = ":".join(keep_keys)
            out_line = toks[:9] + [build_sample_string(keep_keys, dn_out), build_sample_string(keep_keys, rn_out)]
            fout.write("\t".join(out_line) + "\n")
            n_out += 1

    print(f"[done] input_records={n_in} output_records={n_out}", file=sys.stderr)
    print(f"[done] wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
