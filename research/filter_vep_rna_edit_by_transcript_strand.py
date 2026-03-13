#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Quick VEP-only strandedness proxy for RNA-edit-like variants. "
            "Uses VEP transcript STRAND and Uploaded_variation ref/alt to keep only "
            "transcript-canonical ADAR-like A>G or APOBEC3-like C>T substitutions."
        )
    )
    ap.add_argument(
        "--input",
        action="append",
        default=[],
        help="PATIENT=/path/to/patient_vep.vep or just /path/to/patient_vep.vep",
    )
    ap.add_argument("--vep-dir", default="", help="Directory with per-patient VEP files named {patient}_vep.vep")
    ap.add_argument("--patient", action="append", default=[], help="Optional patient id(s) when using --vep-dir")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument(
        "--cohort-prefix",
        default="cohort",
        help="Prefix for merged cohort outputs when multiple inputs are provided",
    )
    return ap.parse_args()


def infer_patient(value):
    if "=" in value:
        patient, path = value.split("=", 1)
        return patient, path
    path = value
    base = os.path.basename(path)
    patient = re.sub(r"_vep\.vep(\.gz)?$", "", base)
    return patient, path


def reverse_complement(base):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return comp.get(base.upper(), base.upper())


def transcript_space_change(ref, alt, strand):
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) != 1 or len(alt) != 1:
        return None
    if strand == "+":
        return f"{ref}>{alt}"
    if strand == "-":
        return f"{reverse_complement(ref)}>{reverse_complement(alt)}"
    return None


def classify_transcript_change(change):
    if change == "A>G":
        return "ADAR"
    if change == "C>T":
        return "APOBEC3"
    return "OTHER"


def parse_uploaded_variation(uploaded):
    text = str(uploaded).strip()
    m = re.match(r"(.+?)_(\d+)_([^/]+)/([^/]+)$", text)
    if not m:
        return None
    chrom, pos, ref, alt = m.groups()
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
    }


def parse_symbol(extra):
    m = re.search(r"(?:^|;)SYMBOL=([^;]+)", str(extra))
    return m.group(1) if m else ""


def parse_strand(extra):
    m = re.search(r"(?:^|;)STRAND=([+-]?[0-9]+)", str(extra))
    if not m:
        return None
    val = m.group(1)
    if val == "1":
        return "+"
    if val == "-1":
        return "-"
    return None


def load_vep(path):
    header = None
    headers_raw = []
    rows = []
    meta = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                meta.append(line)
                continue
            if line.startswith("#"):
                headers_raw = line.rstrip("\n").split("\t")
                header = [c.lstrip("#") for c in headers_raw]
                continue
            cols = line.rstrip("\n").split("\t")
            if not cols or not header:
                continue
            row = {header[i]: cols[i] if i < len(cols) else "" for i in range(len(header))}
            row["_raw_cols"] = cols
            rows.append(row)
    return meta, headers_raw, rows


def evaluate_rows(rows):
    by_uploaded = defaultdict(list)
    for row in rows:
        uploaded = str(row.get("Uploaded_variation", "")).strip()
        if uploaded:
            by_uploaded[uploaded].append(row)

    keep_uploaded = set()
    summary_rows = []

    for uploaded, group in sorted(by_uploaded.items()):
        parsed = parse_uploaded_variation(uploaded)
        if not parsed:
            summary_rows.append(
                {
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
                    "keep": "0",
                    "reason": "UNPARSEABLE_UPLOADED_VARIATION",
                }
            )
            continue

        strand_set = set()
        gene_symbols = set()
        consequences = set()
        feature_types = set()
        locations = set()
        for row in group:
            strand = parse_strand(row.get("Extra", ""))
            if strand:
                strand_set.add(strand)
            symbol = parse_symbol(row.get("Extra", ""))
            if symbol:
                gene_symbols.add(symbol)
            if row.get("Consequence"):
                consequences.add(row["Consequence"])
            if row.get("Feature_type"):
                feature_types.add(row["Feature_type"])
            if row.get("Location"):
                locations.add(row["Location"])

        if not strand_set:
            keep = False
            reason = "MISSING_TRANSCRIPT_STRAND"
            tx_change = "NA"
            edit_class = "NA"
        elif len(strand_set) > 1:
            keep = False
            reason = "AMBIGUOUS_TRANSCRIPT_STRAND"
            tx_change = "NA"
            edit_class = "NA"
        else:
            strand = next(iter(strand_set))
            tx_change = transcript_space_change(parsed["ref"], parsed["alt"], strand)
            if tx_change is None:
                keep = False
                reason = "NON_SNV"
                edit_class = "NA"
            else:
                edit_class = classify_transcript_change(tx_change)
                if edit_class in ("ADAR", "APOBEC3"):
                    keep = True
                    reason = "PASS"
                    keep_uploaded.add(uploaded)
                else:
                    keep = False
                    reason = "NOT_TRANSCRIPT_CANONICAL_RNA_EDIT"

        summary_rows.append(
            {
                "uploaded_variation": uploaded,
                "chrom": parsed["chrom"],
                "pos": parsed["pos"],
                "ref": parsed["ref"],
                "alt": parsed["alt"],
                "location": "|".join(sorted(locations)) if locations else "NA",
                "gene_symbol": "|".join(sorted(gene_symbols)) if gene_symbols else "NA",
                "consequence": "|".join(sorted(consequences)) if consequences else "NA",
                "feature_type": "|".join(sorted(feature_types)) if feature_types else "NA",
                "transcript_strand_set": "|".join(sorted(strand_set)) if strand_set else "NA",
                "transcript_space_change": tx_change,
                "edit_class": edit_class,
                "keep": "1" if keep else "0",
                "reason": reason,
            }
        )

    return keep_uploaded, summary_rows


def write_vep(path, meta, headers_raw, rows, keep_uploaded):
    with open(path, "w") as out:
        for line in meta:
            out.write(line)
        if headers_raw:
            out.write("\t".join(headers_raw) + "\n")
        for row in rows:
            if row.get("Uploaded_variation", "") in keep_uploaded:
                out.write("\t".join(row["_raw_cols"]) + "\n")


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_list(path, values):
    with open(path, "w") as out:
        for value in sorted(set(values)):
            out.write(f"{value}\n")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    inputs = list(args.input)
    if args.vep_dir:
        if args.patient:
            for patient in args.patient:
                inputs.append(f"{patient}={os.path.join(args.vep_dir, f'{patient}_vep.vep')}")
        else:
            for name in sorted(os.listdir(args.vep_dir)):
                if re.search(r"_vep\.vep(\.gz)?$", name):
                    patient = re.sub(r"_vep\.vep(\.gz)?$", "", name)
                    inputs.append(f"{patient}={os.path.join(args.vep_dir, name)}")
    if not inputs:
        raise SystemExit("ERROR: provide --input and/or --vep-dir")

    summary_fieldnames = [
        "uploaded_variation",
        "chrom",
        "pos",
        "ref",
        "alt",
        "location",
        "gene_symbol",
        "consequence",
        "feature_type",
        "transcript_strand_set",
        "transcript_space_change",
        "edit_class",
        "keep",
        "reason",
    ]

    cohort_summary = []
    cohort_keep = []
    cohort_blacklist = []

    for item in inputs:
        patient, path = infer_patient(item)
        if not os.path.exists(path):
            print(f"[warn] missing VEP for {patient}: {path}")
            continue

        meta, headers_raw, rows = load_vep(path)
        keep_uploaded, summary_rows = evaluate_rows(rows)
        patient_summary = []
        for row in summary_rows:
            row2 = dict(row)
            row2["SAMPLE"] = patient
            patient_summary.append(row2)

        patient_prefix = os.path.join(args.outdir, patient)
        filtered_vep = f"{patient_prefix}.vep.transcript_strand_filtered.vep"
        summary_tsv = f"{patient_prefix}.vep.transcript_strand_summary.tsv"
        keep_txt = f"{patient_prefix}.vep.transcript_strand_keep.uploaded_variation.txt"
        blacklist_txt = f"{patient_prefix}.vep.transcript_strand_blacklist.uploaded_variation.txt"

        write_vep(filtered_vep, meta, headers_raw, rows, keep_uploaded)
        write_tsv(summary_tsv, patient_summary, ["SAMPLE"] + summary_fieldnames)
        write_list(keep_txt, keep_uploaded)
        write_list(blacklist_txt, [r["uploaded_variation"] for r in summary_rows if r["keep"] != "1"])

        print(f"[done] {patient}: kept={len(keep_uploaded)} filtered_vep={filtered_vep}")

        cohort_summary.extend(patient_summary)
        cohort_keep.extend([f"{patient}\t{x}" for x in keep_uploaded])
        cohort_blacklist.extend([f"{patient}\t{r['uploaded_variation']}" for r in summary_rows if r["keep"] != "1"])

    if len(inputs) > 1 and cohort_summary:
        cohort_summary_path = os.path.join(args.outdir, f"{args.cohort_prefix}.vep.transcript_strand_summary.tsv")
        cohort_keep_path = os.path.join(args.outdir, f"{args.cohort_prefix}.vep.transcript_strand_keep.tsv")
        cohort_blacklist_path = os.path.join(args.outdir, f"{args.cohort_prefix}.vep.transcript_strand_blacklist.tsv")
        write_tsv(cohort_summary_path, cohort_summary, ["SAMPLE"] + summary_fieldnames)
        with open(cohort_keep_path, "w") as out:
            out.write("SAMPLE\tuploaded_variation\n")
            for line in sorted(set(cohort_keep)):
                out.write(f"{line}\n")
        with open(cohort_blacklist_path, "w") as out:
            out.write("SAMPLE\tuploaded_variation\n")
            for line in sorted(set(cohort_blacklist)):
                out.write(f"{line}\n")
        print(f"[done] cohort summary -> {cohort_summary_path}")


if __name__ == "__main__":
    main()
