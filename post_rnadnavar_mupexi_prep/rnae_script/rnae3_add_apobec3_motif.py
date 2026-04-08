#!/usr/bin/env python3
"""
Annotate APOBEC3 motif matches into INFO/KNOWN_RNAEDIT_DB.

The motif is interpreted in transcript space as CAT[C>T]G.
Using only the reference genome:
  - genomic C>T matches if the + strand context is CATCG
  - genomic G>A matches if the reverse-complement equivalent + strand context is CGATG

The script appends APOBEC3_MOTIF into KNOWN_RNAEDIT_DB for APOBEC3 variants
detected by INFO/APOBEC3_SIG or INFO/EDIT_SIG containing APOBEC3/APOBEC.
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
from typing import Dict, Iterable, Optional, TextIO


def open_in(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def open_out(path: str):
    if path == "-":
        return sys.stdout
    if path.endswith(".gz"):
        try:
            import pysam  # type: ignore
        except ImportError as exc:
            raise SystemExit(
                "ERROR: pysam is required to write bgzipped output. "
                "Load a Python environment with pysam."
            ) from exc
        return io.TextIOWrapper(pysam.BGZFile(path, "w"), encoding="utf-8")
    return open(path, "w")


def build_contig_lookup(references: Iterable[str]) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for ref in references:
        lookup.setdefault(ref, ref)
        if ref.lower().startswith("chr"):
            lookup.setdefault(ref[3:], ref)
        else:
            lookup.setdefault(f"chr{ref}", ref)
    return lookup


def parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if info in ("", "."):
        return out
    for token in info.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
        else:
            out[token] = "1"
    return out


def build_info_string(original_info: str, key: str, label: str) -> str:
    if original_info in ("", "."):
        return f"{key}={label}"

    fields = original_info.split(";")
    found = False
    updated = []
    for field in fields:
        if field.startswith(f"{key}="):
            found = True
            values = [v for v in field.split("=", 1)[1].split(",") if v not in ("", ".")]
            if label not in values:
                values.append(label)
            updated.append(f"{key}={','.join(values)}")
        else:
            updated.append(field)
    if not found:
        updated.append(f"{key}={label}")
    return ";".join(updated)


def is_apobec3_variant(info_map: Dict[str, str]) -> bool:
    if "APOBEC3_SIG" in info_map:
        return True
    edit_sig = str(info_map.get("EDIT_SIG", "")).upper()
    return "APOBEC3" in edit_sig or "APOBEC" in edit_sig


def motif_match(ref: str, alt: str, context_getter, chrom: str, pos: int) -> bool:
    ref = ref.upper()
    alt = alt.upper()
    if ref == "C" and alt == "T":
        seq = context_getter(chrom, pos - 3, pos + 1)
        return seq == "CATCG"
    if ref == "G" and alt == "A":
        seq = context_getter(chrom, pos - 1, pos + 3)
        return seq == "CGATG"
    return False


def main():
    ap = argparse.ArgumentParser(description="Append APOBEC3_MOTIF to INFO/KNOWN_RNAEDIT_DB using reference-genome motif context.")
    ap.add_argument("-i", "--input", required=True, help="Input VCF/VCF.GZ")
    ap.add_argument("-o", "--output", required=True, help="Output VCF/VCF.GZ")
    ap.add_argument("-f", "--fasta", required=True, help="Reference FASTA used for motif lookup")
    ap.add_argument("--info-field", default="KNOWN_RNAEDIT_DB", help="INFO field to append the motif label to")
    ap.add_argument("--label", default="APOBEC3_MOTIF", help="Label value to append")
    args = ap.parse_args()

    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "ERROR: pysam is required for APOBEC3 motif annotation. "
            "Load a Python environment with pysam."
        ) from exc

    fasta = pysam.FastaFile(args.fasta)
    contig_lookup = build_contig_lookup(fasta.references)

    def fetch_context(chrom: str, start_1based: int, end_1based: int) -> Optional[str]:
        actual = contig_lookup.get(chrom)
        if actual is None:
            return None
        if start_1based < 1:
            return None
        try:
            seq = fasta.fetch(actual, start_1based - 1, end_1based)
        except Exception:
            return None
        seq = seq.upper()
        expected_len = end_1based - start_1based + 1
        if len(seq) != expected_len:
            return None
        return seq

    hdr = (
        f'##INFO=<ID={args.info_field},Number=.,Type=String,'
        'Description="Known RNA-editing DB hit(s) and motif-based label(s); merged sources keyed by CHROM,POS,REF,ALT">'
    )

    total = 0
    apobec_total = 0
    motif_hits = 0
    newly_tagged = 0
    header_written = False

    fin = open_in(args.input)
    fout = open_out(args.output)
    try:
        for line in fin:
            if line.startswith(f"##INFO=<ID={args.info_field},"):
                if not header_written:
                    fout.write(hdr + "\n")
                    header_written = True
                continue
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                if not header_written:
                    fout.write(hdr + "\n")
                    header_written = True
                fout.write(line)
                continue
            if not line or line.startswith("#"):
                fout.write(line)
                continue

            total += 1
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 8:
                fout.write(line)
                continue

            chrom = toks[0]
            pos = int(toks[1])
            ref = toks[3]
            alt = toks[4]
            if len(ref) != 1 or len(alt) != 1:
                fout.write(line)
                continue

            info = toks[7]
            info_map = parse_info(info)
            if is_apobec3_variant(info_map):
                apobec_total += 1
                if motif_match(ref, alt, fetch_context, chrom, pos):
                    motif_hits += 1
                    had_label = args.info_field in info_map and args.label in str(info_map.get(args.info_field, "")).split(",")
                    toks[7] = build_info_string(info, args.info_field, args.label)
                    if not had_label:
                        newly_tagged += 1

            fout.write("\t".join(toks) + "\n")
    finally:
        fin.close()
        if fout is not sys.stdout:
            fout.close()
        fasta.close()

    print(f"[done] wrote {args.output}", file=sys.stderr)
    print(f"records_total_in\t{total}", file=sys.stderr)
    print(f"records_apobec3\t{apobec_total}", file=sys.stderr)
    print(f"records_apobec3_motif_hits\t{motif_hits}", file=sys.stderr)
    print(f"records_newly_tagged\t{newly_tagged}", file=sys.stderr)


if __name__ == "__main__":
    main()
