#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, TextIO, Tuple


PREFERRED_METHOD_NOTE = (
    "Post hoc VEP-only strandedness filter applied at rna7.1. "
    "Preferred method is rna5.1_FilterByStrandness.sh before merge/phasing."
)


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom(chrom: str) -> str:
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def parse_uploaded_variation(uploaded: str) -> Optional[Tuple[str, int, str, str]]:
    match = re.match(r"(.+?)_(\d+)_([^/]+)/([^/]+)$", str(uploaded).strip())
    if not match:
        return None
    chrom, pos_s, ref, alt = match.groups()
    try:
        pos = int(pos_s)
    except ValueError:
        return None
    return normalize_chrom(chrom), pos, ref, alt


def reverse_complement(base: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return comp.get(base.upper(), base.upper())


def transcript_space_change(ref: str, alt: str, strand: str) -> Optional[str]:
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) != 1 or len(alt) != 1:
        return None
    if strand == "+":
        return f"{ref}>{alt}"
    if strand == "-":
        return f"{reverse_complement(ref)}>{reverse_complement(alt)}"
    return None


def classify_transcript_change(change: Optional[str]) -> str:
    if change == "A>G":
        return "ADAR"
    if change == "C>T":
        return "APOBEC3"
    return "OTHER"


def parse_symbol(extra: str) -> str:
    match = re.search(r"(?:^|;)SYMBOL=([^;]+)", str(extra))
    return match.group(1) if match else ""


def parse_strand(extra: str) -> Optional[str]:
    match = re.search(r"(?:^|;)STRAND=([+-]?[0-9]+)", str(extra))
    if not match:
        return None
    value = match.group(1)
    if value == "1":
        return "+"
    if value == "-1":
        return "-"
    return None


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


def join_values(values: Iterable[str]) -> str:
    clean = sorted({str(v) for v in values if str(v) not in ("", ".", "NA")})
    return "|".join(clean) if clean else "NA"


def rec_strand_status(decision: Dict[str, str]) -> str:
    return str(decision.get("strand_status", "NA"))


def fallback_uploaded_variation(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{normalize_chrom(chrom)}_{pos}_{ref}/{alt}"


def load_vep_rows(path: str) -> List[Dict[str, str]]:
    header: Optional[List[str]] = None
    rows: List[Dict[str, str]] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = [c.lstrip("#") for c in line.rstrip("\n").split("\t")]
                continue
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if header:
                row = {header[i]: cols[i] if i < len(cols) else "" for i in range(len(header))}
            else:
                default_header = [
                    "Uploaded_variation",
                    "Location",
                    "Allele",
                    "Gene",
                    "Feature",
                    "Feature_type",
                    "Consequence",
                    "cDNA_position",
                    "CDS_position",
                    "Protein_position",
                    "Amino_acids",
                    "Codons",
                    "Existing_variation",
                    "Extra",
                ]
                row = {default_header[i]: cols[i] if i < len(cols) else "" for i in range(len(default_header))}
            rows.append(row)
    return rows


def summarize_vep(rows: Sequence[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    by_uploaded: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        uploaded = str(row.get("Uploaded_variation", "")).strip()
        if uploaded:
            by_uploaded[uploaded].append(row)

    summary: Dict[str, Dict[str, str]] = {}
    for uploaded, group in sorted(by_uploaded.items()):
        parsed = parse_uploaded_variation(uploaded)
        if not parsed:
            summary[uploaded] = {
                "uploaded_variation": uploaded,
                "chrom": "NA",
                "pos": "NA",
                "ref": "NA",
                "alt": "NA",
                "location": "NA",
                "gene_symbol": "NA",
                "consequence": "NA",
                "feature_type": "NA",
                "transcript_strand_set": "NA",
                "transcript_space_change": "NA",
                "edit_class": "NA",
                "strand_status": "UNPARSEABLE_UPLOADED_VARIATION",
            }
            continue

        chrom, pos, ref, alt = parsed
        strand_set: Set[str] = set()
        gene_symbols: Set[str] = set()
        consequences: Set[str] = set()
        feature_types: Set[str] = set()
        locations: Set[str] = set()
        for row in group:
            strand = parse_strand(row.get("Extra", ""))
            if strand:
                strand_set.add(strand)
            symbol = parse_symbol(row.get("Extra", ""))
            if symbol:
                gene_symbols.add(symbol)
            value = str(row.get("Consequence", "")).strip()
            if value:
                consequences.add(value)
            value = str(row.get("Feature_type", "")).strip()
            if value:
                feature_types.add(value)
            value = str(row.get("Location", "")).strip()
            if value:
                locations.add(value)

        if not strand_set:
            strand_status = "MISSING_TRANSCRIPT_STRAND"
            tx_change = "NA"
            edit_class = "NA"
        elif len(strand_set) > 1:
            strand_status = "AMBIGUOUS_TRANSCRIPT_STRAND"
            tx_change = "NA"
            edit_class = "NA"
        else:
            strand = next(iter(strand_set))
            tx_change = transcript_space_change(ref, alt, strand)
            strand_status = "UNIQUE_TRANSCRIPT_STRAND"
            if tx_change is None:
                strand_status = "NON_SNV"
                edit_class = "NA"
            else:
                edit_class = classify_transcript_change(tx_change)

        summary[uploaded] = {
            "uploaded_variation": uploaded,
            "chrom": chrom,
            "pos": str(pos),
            "ref": ref,
            "alt": alt,
            "location": join_values(locations),
            "gene_symbol": join_values(gene_symbols),
            "consequence": join_values(consequences),
            "feature_type": join_values(feature_types),
            "transcript_strand_set": join_values(strand_set),
            "transcript_space_change": tx_change,
            "edit_class": edit_class,
            "strand_status": strand_status,
        }

    return summary


def write_tsv(path: str, rows: Sequence[Dict[str, str]], fieldnames: Sequence[str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_stats(path: str, rows: Sequence[Tuple[str, str]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in rows:
            writer.writerow([key, value])


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Post hoc VEP-only RNA strandedness filter for final phased VCFs. "
            "Preferred method is rna5.1 before merge/phasing."
        )
    )
    ap.add_argument("-i", "--input", required=True, help="Input phased VCF (preferably already split to one ALT per row)")
    ap.add_argument("-o", "--output", required=True, help="Output filtered VCF (plain text; bgzip outside if needed)")
    ap.add_argument("--vep", required=True, help="Per-patient VEP file")
    ap.add_argument("--patient", default="", help="Patient/sample id for sidecar outputs")
    ap.add_argument("--keep-tsv", required=True, help="Output TSV with retained RNA_EDIT variants")
    ap.add_argument("--blacklist-tsv", required=True, help="Output TSV with filtered RNA_EDIT variants")
    ap.add_argument("--stats", required=True, help="Output summary TSV")
    args = ap.parse_args()

    vep_summary = summarize_vep(load_vep_rows(args.vep))

    keep_rows: List[Dict[str, str]] = []
    blacklist_rows: List[Dict[str, str]] = []
    counts: Counter[str] = Counter()

    def flush_position_group(records: Sequence[Dict[str, object]], fout: TextIO) -> None:
        if not records:
            return

        counts["position_groups_total"] += 1
        rna_edit_indices = [idx for idx, rec in enumerate(records) if rec["is_rna_edit"]]
        if not rna_edit_indices:
            for rec in records:
                counts["total_records"] += 1
                counts["non_rna_edit_records_retained"] += 1
                fout.write(rec["line"])
            return

        counts["position_groups_with_rna_edit"] += 1
        for _ in rna_edit_indices:
            counts["rna_edit_records_total"] += 1

        filtered_indices: Set[int] = set()
        keep_reason_by_index: Dict[int, str] = {}
        filter_reason_by_index: Dict[int, str] = {}

        if len(rna_edit_indices) == 1:
            keep_reason_by_index[rna_edit_indices[0]] = "SINGLETON_KEEP"
        else:
            counts["position_groups_with_multiple_rna_edit"] += 1
            unique_indices = [
                idx
                for idx in rna_edit_indices
                if rec_strand_status(records[idx]["decision"]) == "UNIQUE_TRANSCRIPT_STRAND"
            ]
            non_unique_indices = [idx for idx in rna_edit_indices if idx not in unique_indices]

            if len(unique_indices) == 1 and non_unique_indices:
                counts["resolved_position_conflict_groups"] += 1
                keep_reason_by_index[unique_indices[0]] = "UNIQUE_STRAND_KEEP_IN_POSITION_CONFLICT"
                for idx in non_unique_indices:
                    filtered_indices.add(idx)
                    filter_reason_by_index[idx] = "POSITION_CONFLICT_NONUNIQUE_STRAND"
            else:
                counts["unresolved_position_conflict_groups"] += 1
                for idx in rna_edit_indices:
                    keep_reason_by_index[idx] = "UNRESOLVED_POSITION_CONFLICT_KEEP_ALL"

        for idx, rec in enumerate(records):
            counts["total_records"] += 1
            if not rec["is_rna_edit"]:
                counts["non_rna_edit_records_retained"] += 1
                fout.write(rec["line"])
                continue

            decision = rec["decision"]
            row = {
                "SAMPLE": args.patient or "NA",
                "uploaded_variation": rec["uploaded"],
                "chrom": rec["chrom"],
                "pos": str(rec["pos"]),
                "ref": rec["ref"],
                "alt": rec["alt"],
                "source_set": join_values(rec["source_tokens"]),
                "filter": rec["filter"],
                "edit_class": decision["edit_class"],
                "transcript_strand_set": decision["transcript_strand_set"],
                "transcript_space_change": decision["transcript_space_change"],
                "strand_status": rec_strand_status(decision),
                "position_group_size": str(len(rna_edit_indices)),
                "reason": keep_reason_by_index.get(idx, filter_reason_by_index.get(idx, "NA")),
            }

            if idx in filtered_indices:
                counts["rna_edit_records_filtered"] += 1
                counts[f"reason::{row['reason']}"] += 1
                blacklist_rows.append(row)
            else:
                counts["rna_edit_records_retained"] += 1
                counts[f"reason::{row['reason']}"] += 1
                keep_rows.append(row)
                fout.write(rec["line"])

    def fallback_decision(uploaded: str, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, str]:
        return {
            "uploaded_variation": uploaded,
            "chrom": chrom,
            "pos": str(pos),
            "ref": ref,
            "alt": alt,
            "location": f"{chrom}:{pos}",
            "gene_symbol": "NA",
            "consequence": "NA",
            "feature_type": "NA",
            "transcript_strand_set": "NA",
            "transcript_space_change": "NA",
            "edit_class": "NA",
            "strand_status": "MISSING_VEP_ANNOTATION",
        }

    def build_record(line: str) -> Dict[str, object]:
        cols = line.rstrip("\n").split("\t")
        chrom = normalize_chrom(cols[0])
        pos = int(cols[1])
        ref = cols[3]
        alt = cols[4]
        if "," in alt:
            raise SystemExit(
                "ERROR: rna7.1 input contains multiallelic records. "
                "Split the VCF first (bcftools norm -m -any)."
            )
        info_map = parse_info(cols[7])
        source_tokens = {token.strip() for token in str(info_map.get("SOURCE_SET", "")).split(",") if token.strip()}
        uploaded = fallback_uploaded_variation(chrom, pos, ref, alt)
        decision = vep_summary.get(uploaded)
        if decision is None:
            decision = fallback_decision(uploaded, chrom, pos, ref, alt)
        return {
            "line": line,
            "cols": cols,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "filter": cols[6],
            "uploaded": uploaded,
            "source_tokens": source_tokens,
            "decision": decision,
            "is_rna_edit": "RNA_EDIT" in source_tokens,
        }

    with open_text(args.input) as fin, open(args.output, "w") as fout:
        inserted_note = False
        current_key: Optional[Tuple[str, int]] = None
        current_records: List[Dict[str, object]] = []
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                if not inserted_note:
                    fout.write(f"##rnadnavar_posthoc_strand_filter={PREFERRED_METHOD_NOTE}\n")
                    inserted_note = True
                fout.write(line)
                continue
            if not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            key = (normalize_chrom(cols[0]), int(cols[1]))
            if current_key is not None and key != current_key:
                flush_position_group(current_records, fout)
                current_records = []
            current_key = key
            current_records.append(build_record(line))

        flush_position_group(current_records, fout)

    fieldnames = [
        "SAMPLE",
        "uploaded_variation",
        "chrom",
        "pos",
        "ref",
        "alt",
        "source_set",
        "filter",
        "edit_class",
        "transcript_strand_set",
        "transcript_space_change",
        "strand_status",
        "position_group_size",
        "reason",
    ]
    write_tsv(args.keep_tsv, keep_rows, fieldnames)
    write_tsv(args.blacklist_tsv, blacklist_rows, fieldnames)

    stats_rows: List[Tuple[str, str]] = [
        ("patient", args.patient or "NA"),
        ("input_vcf", os.path.abspath(args.input)),
        ("output_vcf", os.path.abspath(args.output)),
        ("vep", os.path.abspath(args.vep)),
        ("note", PREFERRED_METHOD_NOTE),
        ("total_records", str(counts["total_records"])),
        ("non_rna_edit_records_retained", str(counts["non_rna_edit_records_retained"])),
        ("position_groups_total", str(counts["position_groups_total"])),
        ("position_groups_with_rna_edit", str(counts["position_groups_with_rna_edit"])),
        ("position_groups_with_multiple_rna_edit", str(counts["position_groups_with_multiple_rna_edit"])),
        ("resolved_position_conflict_groups", str(counts["resolved_position_conflict_groups"])),
        ("unresolved_position_conflict_groups", str(counts["unresolved_position_conflict_groups"])),
        ("rna_edit_records_total", str(counts["rna_edit_records_total"])),
        ("rna_edit_records_retained", str(counts["rna_edit_records_retained"])),
        ("rna_edit_records_filtered", str(counts["rna_edit_records_filtered"])),
    ]
    reason_keys = sorted(key for key in counts if key.startswith("reason::"))
    for key in reason_keys:
        stats_rows.append((key.replace("reason::", "reason_count."), str(counts[key])))
    write_stats(args.stats, stats_rows)

    print(
        "[info] rna7.1 posthoc strand filter: "
        f"rna_edit_total={counts['rna_edit_records_total']} "
        f"retained={counts['rna_edit_records_retained']} "
        f"filtered={counts['rna_edit_records_filtered']}"
    )
    print(f"[warn] preferred method remains rna5.1_FilterByStrandness.sh before merge/phasing")


if __name__ == "__main__":
    main()
