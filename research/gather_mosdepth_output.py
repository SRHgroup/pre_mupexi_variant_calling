#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
from typing import List


def sample_base_name(value: str, labels: List[str]) -> str:
    out = value
    for label in labels:
        suffix = f"_{label}"
        if out.endswith(suffix):
            out = out[: -len(suffix)]
    return out


def load_patients(samples_path: str, labels: List[str], only_patient: str) -> List[str]:
    seen = set()
    patients: List[str] = []
    with open(samples_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            sid = re.split(r"[,\t ]+", line)[0]
            if sid.lower() in {"sample", "sample_id", "patient", "patient_id", "#sample"}:
                continue
            patient = sample_base_name(sid, labels)
            if not patient:
                continue
            if only_patient and only_patient not in (sid, patient):
                continue
            if patient in seen:
                continue
            seen.add(patient)
            patients.append(patient)
    return patients


def find_summary_file(mosdepth_dir: str, patient: str, label: str) -> str:
    candidates = [
        os.path.join(mosdepth_dir, f"{patient}_{label}", f"{patient}_{label}.md.mosdepth.summary.txt"),
        os.path.join(mosdepth_dir, f"{patient}_{patient}_{label}", f"{patient}_{patient}_{label}.md.mosdepth.summary.txt"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path

    patterns = [
        os.path.join(mosdepth_dir, f"{patient}_{label}", "*.mosdepth.summary.txt"),
        os.path.join(mosdepth_dir, f"{patient}_{patient}_{label}", "*.mosdepth.summary.txt"),
        os.path.join(mosdepth_dir, f"{patient}_*_{label}", "*.mosdepth.summary.txt"),
        os.path.join(mosdepth_dir, f"{patient}*{label}", "*.mosdepth.summary.txt"),
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    return ""


def extract_total_mean(path: str) -> str:
    if not path or not os.path.isfile(path):
        return "NA"
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if parts[0] == "total":
                if len(parts) < 4:
                    return "NA"
                return parts[3]
    return "NA"


def main() -> None:
    ap = argparse.ArgumentParser(description="Gather mosdepth total mean coverage for DNA normal, DNA tumor, and RNA tumor.")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--mosdepth-dir", required=True)
    ap.add_argument("--dna-normal-label", required=True)
    ap.add_argument("--dna-tumor-label", required=True)
    ap.add_argument("--rna-tumor-label", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--patient", default="")
    args = ap.parse_args()

    labels = [args.dna_normal_label, args.dna_tumor_label, args.rna_tumor_label, "DNA_NORMAL", "DNA_TUMOR", "DNA_TUMOUR", "RNA_TUMOR", "RNA_TUMOUR", "TUMOR", "TUMOUR"]
    patients = load_patients(args.samples, labels, args.patient)
    if not patients:
        raise SystemExit("ERROR: no patients found in samples file")

    rows = []
    for patient in patients:
        dn_file = find_summary_file(args.mosdepth_dir, patient, args.dna_normal_label)
        dt_file = find_summary_file(args.mosdepth_dir, patient, args.dna_tumor_label)
        rt_file = find_summary_file(args.mosdepth_dir, patient, args.rna_tumor_label)
        rows.append(
            {
                "patient_id": patient,
                "dna_normal_total_mean_coverage": extract_total_mean(dn_file),
                "dna_tumor_total_mean_coverage": extract_total_mean(dt_file),
                "rna_tumor_total_mean_coverage": extract_total_mean(rt_file),
                "dna_normal_summary_file": dn_file or "NA",
                "dna_tumor_summary_file": dt_file or "NA",
                "rna_tumor_summary_file": rt_file or "NA",
            }
        )

    rows.sort(key=lambda r: r["patient_id"])
    fieldnames = [
        "patient_id",
        "dna_normal_total_mean_coverage",
        "dna_tumor_total_mean_coverage",
        "rna_tumor_total_mean_coverage",
        "dna_normal_summary_file",
        "dna_tumor_summary_file",
        "rna_tumor_summary_file",
    ]
    os.makedirs(os.path.dirname(args.outfile) or ".", exist_ok=True)
    with open(args.outfile, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"[done] wrote {len(rows)} mosdepth coverage rows -> {args.outfile}")


if __name__ == "__main__":
    main()
