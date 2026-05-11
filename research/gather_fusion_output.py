#!/usr/bin/env python3
import argparse
import csv
import os
from collections import OrderedDict


STANDARD_FIELDS = [
    "patient_id",
    "event_type",
    "event_id",
    "source_set",
    "source_file",
    "fusion_id",
    "gene1",
    "gene2",
    "breakpoint1",
    "breakpoint2",
]


def parse_map(items):
    out = []
    for item in items:
        if "=" not in item:
            raise SystemExit(f"ERROR: expected PATIENT=/path form, got: {item}")
        patient, path = item.split("=", 1)
        out.append((patient, path))
    return out


def load_arriba_rows(path):
    header = None
    rows = []
    with open(path, "r", newline="") as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            row = {header[i]: cols[i] if i < len(cols) else "" for i in range(len(header))}
            rows.append(row)
    return header or [], rows


def build_fusion_id(row):
    gene1 = row.get("#gene1", row.get("gene1", ""))
    gene2 = row.get("gene2", "")
    breakpoint1 = row.get("breakpoint1", "")
    breakpoint2 = row.get("breakpoint2", "")
    if gene1 and gene2 and breakpoint1 and breakpoint2:
        return f"{gene1}|{gene2}_{breakpoint1}|{breakpoint2}"
    return "NA"


def gather(inputs):
    rows = []
    fieldnames = []
    seen = OrderedDict()
    for patient, path in inputs:
        header, file_rows = load_arriba_rows(path)
        if header and not fieldnames:
            fieldnames = STANDARD_FIELDS + [c for c in header if c not in STANDARD_FIELDS]
        for row in file_rows:
            out = dict(row)
            out["patient_id"] = patient
            out["event_type"] = "FUS"
            out["source_set"] = "FUSION"
            out["source_file"] = path
            out["gene1"] = row.get("#gene1", row.get("gene1", ""))
            out["gene2"] = row.get("gene2", "")
            out["breakpoint1"] = row.get("breakpoint1", "")
            out["breakpoint2"] = row.get("breakpoint2", "")
            out["fusion_id"] = build_fusion_id(row)
            out["event_id"] = out["fusion_id"]
            seen[(patient, out["event_id"])] = out
    rows = list(seen.values())
    rows.sort(key=lambda r: (r.get("patient_id", ""), r.get("event_id", "")))
    if not fieldnames:
        fieldnames = STANDARD_FIELDS
    return fieldnames, rows


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    ap = argparse.ArgumentParser(description="Gather Arriba fusion TSVs into a cohort summary with MuPeXI-compatible fusion_id values.")
    ap.add_argument("--fusion-input", action="append", required=True, help="PATIENT=/path/to/fusions_arriba.tsv")
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    inputs = parse_map(args.fusion_input)
    fieldnames, rows = gather(inputs)
    write_tsv(args.outfile, rows, fieldnames)
    print(f"[done] wrote {len(rows)} fusion rows -> {args.outfile}")


if __name__ == "__main__":
    main()
