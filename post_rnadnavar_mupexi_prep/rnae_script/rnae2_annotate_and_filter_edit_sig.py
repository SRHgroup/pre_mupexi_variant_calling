#!/usr/bin/env python3
"""
annotate_and_filter_edit_sig.py

Annotate RNA-editing substitution signatures and REMOVE all records
that do not match the selected signature(s).

Signatures:
  ADAR-like:    A>G or T>C
  APOBEC3-like: C>T or G>A

Output keeps only SNVs by default (REF and ALT alleles length 1).
For multi-allelic records, the record is kept if ANY ALT allele matches.
(Recommended: split first with `bcftools norm -m -both` for allele-precise filtering.)

Usage:
  python3 annotate_and_filter_edit_sig.py -i in.vcf.gz -o out.vcf.gz
  python3 annotate_and_filter_edit_sig.py -i in.vcf.gz -o out.vcf.gz --adar-only
  python3 annotate_and_filter_edit_sig.py -i in.vcf.gz -o out.vcf.gz --include-apobec3
"""

from __future__ import annotations
import argparse, gzip, sys
from typing import TextIO, Tuple

def open_in(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def open_out(path: str) -> TextIO:
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

def is_snv(ref: str, alt: str) -> bool:
    if len(ref) != 1:
        return False
    return all(len(a) == 1 for a in alt.split(","))

def classify(ref: str, alt: str) -> Tuple[bool, bool]:
    """
    Returns (adar, apobec3) for the record.
    For multi-ALT: true if any ALT allele matches.
    """
    adar = False
    apob = False
    for a in alt.split(","):
        if len(ref) == 1 and len(a) == 1:
            if (ref == "A" and a == "G") or (ref == "T" and a == "C"):
                adar = True
            if (ref == "C" and a == "T") or (ref == "G" and a == "A"):
                apob = True
    return adar, apob

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input VCF (.vcf or .vcf.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output VCF (.vcf or .vcf.gz)")
    ap.add_argument("--adar-only", action="store_true",
                    help="Keep only ADAR-like (A>G/T>C). Default keeps ADAR and APOBEC3.")
    ap.add_argument("--include-apobec3", action="store_true",
                    help="Explicitly include APOBEC3-like (C>T/G>A). If --adar-only is set, ignored.")
    ap.add_argument("--keep-non-snv", action="store_true",
                    help="Do not drop non-SNV records (not recommended). Default: drop non-SNV.")
    args = ap.parse_args()

    # decide what to keep
    keep_adar = True
    keep_apob = (not args.adar_only) and (args.include_apobec3 or True)  # default True unless adar-only

    adar_hdr = '##INFO=<ID=ADAR_SIG,Number=0,Type=Flag,Description="Substitution matches ADAR A->I pattern (A>G or T>C)">\n'
    apob_hdr = '##INFO=<ID=APOBEC3_SIG,Number=0,Type=Flag,Description="Substitution matches APOBEC3 C->U pattern (C>T or G>A)">\n'

    seen_adar = False
    seen_apob = False

    n_total = 0
    n_written = 0
    n_drop_nonsnv = 0
    n_drop_nosig = 0
    n_adar = 0
    n_apob = 0
    n_both = 0

    with open_in(args.input) as fin, open_out(args.output) as fout:
        for line in fin:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=ADAR_SIG"):
                    seen_adar = True
                if line.startswith("##INFO=<ID=APOBEC3_SIG"):
                    seen_apob = True
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                if not seen_adar:
                    fout.write(adar_hdr)
                if not seen_apob:
                    fout.write(apob_hdr)
                fout.write(line)
                continue

            if not line or line[0] == "#":
                fout.write(line)
                continue

            n_total += 1
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 8:
                # malformed; drop
                n_drop_nosig += 1
                continue

            ref = toks[3]
            alt = toks[4]

            if (not args.keep_non_snv) and (not is_snv(ref, alt)):
                n_drop_nonsnv += 1
                continue

            adar, apob = classify(ref, alt)

            # apply keep rules
            keep = (keep_adar and adar) or (keep_apob and apob)
            if not keep:
                n_drop_nosig += 1
                continue

            # annotate INFO flags
            flags = []
            if adar:
                flags.append("ADAR_SIG")
            if apob:
                flags.append("APOBEC3_SIG")

            info = toks[7]
            if info == ".":
                info = ""
            if info:
                info = info + ";" + ";".join(flags)
            else:
                info = ";".join(flags)
            toks[7] = info

            if adar and apob:
                n_both += 1
            elif adar:
                n_adar += 1
            elif apob:
                n_apob += 1

            fout.write("\t".join(toks) + "\n")
            n_written += 1

    print(f"[done] wrote {args.output}", file=sys.stderr)
    print(f"records_total_in\t{n_total}", file=sys.stderr)
    print(f"records_written\t{n_written}", file=sys.stderr)
    print(f"dropped_non_snv\t{n_drop_nonsnv}", file=sys.stderr)
    print(f"dropped_no_signature\t{n_drop_nosig}", file=sys.stderr)
    print(f"written_ADAR_SIG_only\t{n_adar}", file=sys.stderr)
    print(f"written_APOBEC3_SIG_only\t{n_apob}", file=sys.stderr)
    print(f"written_both_flags\t{n_both}", file=sys.stderr)

if __name__ == "__main__":
    main()
