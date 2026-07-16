#!/usr/bin/env python3
"""Gather MuPeXI output files matching a user-specified suffix."""

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


PREFIX_FIELDS = ["patient_id", "source_suffix", "source_file", "source_basename", "source_row_number"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Merge MuPeXI end-product files such as _snv.mupexi, _fus.mupexi, "
            "_neospl.mupexi, or _neojunctions.mupexi into one patient-labelled TSV."
        )
    )
    ap.add_argument("--input-dir", required=True, help="MuPeXI output root to scan recursively")
    ap.add_argument("--suffix", required=True, help="Exact file suffix to gather, e.g. _neospl.mupexi")
    ap.add_argument("--outfile", required=True, help="Merged output TSV")
    ap.add_argument("--patient", action="append", default=[], help="Optional patient filter, repeatable")
    return ap.parse_args()


def normalize_suffix(suffix: str) -> str:
    suffix = suffix.strip()
    if not suffix:
        raise SystemExit("ERROR: --suffix cannot be empty")
    if not suffix.endswith(".mupexi"):
        raise SystemExit(f"ERROR: --suffix must end with .mupexi, got: {suffix}")
    return suffix


def patient_from_path(path: Path, suffix: str) -> str:
    name = path.name
    if not name.endswith(suffix):
        return path.stem
    patient = name[: -len(suffix)]
    return patient.rstrip("._-") or patient


def discover_inputs(root: Path, suffix: str, patients: Sequence[str]) -> List[Path]:
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"ERROR: input dir does not exist or is not a directory: {root}")
    patient_set = set(patients)
    files = []
    for path in root.rglob(f"*{suffix}"):
        if not path.is_file():
            continue
        if path.name.startswith("."):
            continue
        patient = patient_from_path(path, suffix)
        if patient_set and patient not in patient_set:
            continue
        files.append(path)
    return sorted(files)


def iter_mupexi_rows(path: Path) -> Iterable[Dict[str, str]]:
    header = None
    with path.open("r", encoding="utf-8", newline="") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            yield {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}


def extend_fields(fieldnames: List[str], row: Dict[str, str]) -> None:
    seen = set(fieldnames)
    for key in row:
        if key not in seen:
            fieldnames.append(key)
            seen.add(key)


def gather(files: Sequence[Path], suffix: str) -> Tuple[List[str], List[Dict[str, str]]]:
    fieldnames = list(PREFIX_FIELDS)
    rows: List[Dict[str, str]] = []
    for path in files:
        patient = patient_from_path(path, suffix)
        for row_number, row in enumerate(iter_mupexi_rows(path), start=1):
            out = dict(row)
            out.update(
                {
                    "patient_id": patient,
                    "source_suffix": suffix,
                    "source_file": str(path),
                    "source_basename": path.name,
                    "source_row_number": str(row_number),
                }
            )
            extend_fields(fieldnames, out)
            rows.append(out)
    return fieldnames, rows


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    suffix = normalize_suffix(args.suffix)
    input_dir = Path(args.input_dir)
    outfile = Path(args.outfile)
    files = discover_inputs(input_dir, suffix, args.patient)
    if not files:
        raise SystemExit(f"ERROR: no files ending with {suffix} found under {input_dir}")
    fieldnames, rows = gather(files, suffix)
    write_tsv(outfile, fieldnames, rows)
    print(f"[done] suffix={suffix} files={len(files)} rows={len(rows)} -> {outfile}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
