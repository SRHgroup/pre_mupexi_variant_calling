#!/usr/bin/env python3
"""Build a compact normal splice-junction reference from Snaptron GTEx files."""

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


JUNCTION_COLUMNS = [
    "snaptron_id",
    "chrom",
    "start",
    "end",
    "length",
    "strand",
    "annotated",
    "left_motif",
    "right_motif",
    "left_annotated",
    "right_annotated",
    "samples",
    "sample_count",
    "coverage_sum",
    "coverage_avg",
    "coverage_median",
    "source_dataset_id",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Convert a Snaptron GTEx folder containing junctions.bgz, samples.tsv, "
            "and samples.fields.tsv into the compact normal junction TSV used by spl3.5."
        )
    )
    ap.add_argument("--snaptron-dir", required=True, help="Folder containing junctions.bgz, samples.tsv, samples.fields.tsv")
    ap.add_argument("--out", required=True, help="Output compact reference TSV or TSV.GZ")
    ap.add_argument(
        "--coordinate-mode",
        choices=("star-intron", "boundary"),
        default="star-intron",
        help=(
            "How to convert Snaptron start/end. star-intron writes left_boundary=start-1 "
            "and right_boundary=end to match this pipeline's STAR-derived coordinates. "
            "boundary writes left_boundary=start and right_boundary=end."
        ),
    )
    ap.add_argument(
        "--total-samples",
        type=float,
        default=0.0,
        help="Override normal sample denominator. Default: count rows in samples.tsv.",
    )
    ap.add_argument("--min-sample-count", type=float, default=1.0, help="Minimum Snaptron sample_count to retain")
    ap.add_argument("--min-total-reads", type=float, default=1.0, help="Minimum coverage_sum to retain")
    ap.add_argument("--min-prevalence", type=float, default=0.0, help="Minimum normal prevalence to retain")
    ap.add_argument("--canonical-only", action="store_true", help="Keep only chr1-22,X,Y,MT/M")
    ap.add_argument("--drop-unknown-strand", action="store_true", help="Drop junctions with strand ? or missing")
    ap.add_argument(
        "--sample-filter-column",
        action="append",
        default=[],
        metavar="COLUMN=VALUE",
        help="Optional samples.tsv filter used when counting denominator, repeatable. Exact string match.",
    )
    ap.add_argument(
        "--write-summary",
        default="",
        help="Optional summary TSV. Default: OUT with .summary.tsv appended before .gz if OUT ends in .gz.",
    )
    return ap.parse_args()


def open_input(path: Path):
    if path.suffix in {".gz", ".bgz"}:
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8", newline="")


def open_output(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return path.open("w", encoding="utf-8", newline="")


def parse_float(value: str, default: float = 0.0) -> float:
    try:
        return float(str(value).strip())
    except ValueError:
        return default


def chrom_key(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def normalized_chrom(chrom: str) -> str:
    key = chrom_key(chrom)
    if key == "MT":
        return "chrM"
    if key in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        return f"chr{key}"
    return chrom


def is_canonical_chrom(chrom: str) -> bool:
    key = chrom_key(chrom)
    return key in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def parse_sample_filters(raw_filters: Sequence[str]) -> List[Tuple[str, str]]:
    parsed = []
    for raw in raw_filters:
        if "=" not in raw:
            raise SystemExit(f"ERROR: --sample-filter-column must be COLUMN=VALUE, got: {raw}")
        column, value = raw.split("=", 1)
        if not column:
            raise SystemExit(f"ERROR: empty sample filter column in: {raw}")
        parsed.append((column, value))
    return parsed


def load_field_names(fields_path: Path) -> List[str]:
    if not fields_path.exists():
        return []
    names = []
    with fields_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and "field" in reader.fieldnames:
            for row in reader:
                field = row.get("field", "")
                if field:
                    names.append(field)
    return names


def sample_row_matches(row: Dict[str, str], filters: Sequence[Tuple[str, str]]) -> bool:
    for column, value in filters:
        if row.get(column, "") != value:
            return False
    return True


def count_samples(samples_path: Path, fields_path: Path, filters: Sequence[Tuple[str, str]]) -> Tuple[int, Counter]:
    if not samples_path.exists():
        print(f"[snaptron-ref][warn] samples.tsv not found: {samples_path}; total samples denominator will be 0", file=sys.stderr)
        return 0, Counter()

    fields_from_sidecar = load_field_names(fields_path)
    tissue_counts: Counter = Counter()
    count = 0
    with samples_path.open("r", encoding="utf-8", newline="") as fh:
        first_line = fh.readline()
        if not first_line:
            return 0, tissue_counts
        first_fields = first_line.rstrip("\n").split("\t")
        has_header = "rail_id" in first_fields
        fh.seek(0)
        if has_header:
            reader = csv.DictReader(fh, delimiter="\t")
        else:
            if not fields_from_sidecar:
                raise SystemExit("ERROR: samples.tsv has no header and samples.fields.tsv did not provide field names")
            reader = csv.DictReader(fh, delimiter="\t", fieldnames=fields_from_sidecar)
        for row in reader:
            if not sample_row_matches(row, filters):
                continue
            count += 1
            tissue = row.get("SMTSD") or row.get("SMTSC") or row.get("SMTS") or "NA"
            tissue_counts[tissue] += 1
    return count, tissue_counts


def split_junction_line(line: str, line_number: int) -> Optional[Dict[str, str]]:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < len(JUNCTION_COLUMNS):
        print(f"[snaptron-ref][warn] skipping short junction row line={line_number} columns={len(parts)}", file=sys.stderr)
        return None
    return dict(zip(JUNCTION_COLUMNS, parts[: len(JUNCTION_COLUMNS)]))


def output_summary_path(out_path: Path, explicit: str) -> Path:
    if explicit:
        return Path(explicit)
    if out_path.name.endswith(".gz"):
        return out_path.with_name(f"{out_path.name[:-3]}.summary.tsv")
    return out_path.with_name(f"{out_path.name}.summary.tsv")


def write_summary(path: Path, stats: Counter, total_samples: float, tissue_counts: Counter, args: argparse.Namespace) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["snaptron_dir", args.snaptron_dir])
        writer.writerow(["coordinate_mode", args.coordinate_mode])
        writer.writerow(["total_normal_samples", f"{total_samples:g}"])
        for key in sorted(stats):
            writer.writerow([key, stats[key]])
        writer.writerow(["tissue_count_distinct", len(tissue_counts)])
        for tissue, count in tissue_counts.most_common():
            writer.writerow([f"tissue_samples:{tissue}", count])


def build_reference(args: argparse.Namespace) -> Counter:
    snaptron_dir = Path(args.snaptron_dir)
    junctions_path = snaptron_dir / "junctions.bgz"
    samples_path = snaptron_dir / "samples.tsv"
    fields_path = snaptron_dir / "samples.fields.tsv"
    out_path = Path(args.out)
    summary_path = output_summary_path(out_path, args.write_summary)
    filters = parse_sample_filters(args.sample_filter_column)

    if not junctions_path.exists():
        raise SystemExit(f"ERROR: missing {junctions_path}")
    if not fields_path.exists():
        print(f"[snaptron-ref][warn] missing {fields_path}; only needed if samples.tsv has no header", file=sys.stderr)

    counted_samples, tissue_counts = count_samples(samples_path, fields_path, filters)
    total_samples = args.total_samples if args.total_samples > 0 else float(counted_samples)
    if total_samples <= 0:
        print("[snaptron-ref][warn] total sample denominator is 0; normal_prevalence will be NA", file=sys.stderr)

    stats: Counter = Counter()
    fieldnames = [
        "chrom",
        "left_boundary",
        "right_boundary",
        "strand",
        "normal_sample_count",
        "normal_total_reads",
        "normal_prevalence",
        "source",
        "snaptron_id",
        "snaptron_start",
        "snaptron_end",
        "snaptron_length",
        "motif",
        "snaptron_annotated",
        "snaptron_left_annotated",
        "snaptron_right_annotated",
        "normal_avg_reads",
        "normal_median_reads",
        "source_dataset_id",
    ]

    print(f"[snaptron-ref] junctions: {junctions_path}", file=sys.stderr)
    print(f"[snaptron-ref] samples: {samples_path}", file=sys.stderr)
    print(f"[snaptron-ref] total_normal_samples={total_samples:g}", file=sys.stderr)
    print(f"[snaptron-ref] output: {out_path}", file=sys.stderr)

    with open_input(junctions_path) as fh, open_output(out_path) as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line_number, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            stats["input_rows"] += 1
            row = split_junction_line(line, line_number)
            if row is None:
                stats["skipped_short_row"] += 1
                continue

            chrom = row["chrom"]
            if args.canonical_only and not is_canonical_chrom(chrom):
                stats["skipped_noncanonical_chrom"] += 1
                continue
            strand = row["strand"] or "."
            if strand in {"?", "NA"}:
                if args.drop_unknown_strand:
                    stats["skipped_unknown_strand"] += 1
                    continue
                strand = "."

            start = int(row["start"])
            end = int(row["end"])
            left_boundary = start - 1 if args.coordinate_mode == "star-intron" else start
            right_boundary = end
            sample_count = parse_float(row["sample_count"])
            total_reads = parse_float(row["coverage_sum"])
            prevalence = sample_count / total_samples if total_samples > 0 else None
            if sample_count < args.min_sample_count:
                stats["skipped_low_sample_count"] += 1
                continue
            if total_reads < args.min_total_reads:
                stats["skipped_low_total_reads"] += 1
                continue
            if prevalence is not None and prevalence < args.min_prevalence:
                stats["skipped_low_prevalence"] += 1
                continue

            writer.writerow(
                {
                    "chrom": normalized_chrom(chrom),
                    "left_boundary": left_boundary,
                    "right_boundary": right_boundary,
                    "strand": strand,
                    "normal_sample_count": f"{sample_count:g}",
                    "normal_total_reads": f"{total_reads:g}",
                    "normal_prevalence": "NA" if prevalence is None else f"{prevalence:.8g}",
                    "source": "Snaptron_GTEx",
                    "snaptron_id": row["snaptron_id"],
                    "snaptron_start": start,
                    "snaptron_end": end,
                    "snaptron_length": row["length"],
                    "motif": f"{row['left_motif']}/{row['right_motif']}",
                    "snaptron_annotated": row["annotated"],
                    "snaptron_left_annotated": row["left_annotated"],
                    "snaptron_right_annotated": row["right_annotated"],
                    "normal_avg_reads": row["coverage_avg"],
                    "normal_median_reads": row["coverage_median"],
                    "source_dataset_id": row["source_dataset_id"],
                }
            )
            stats["output_rows"] += 1
            if stats["input_rows"] % 1_000_000 == 0:
                print(
                    f"[snaptron-ref] processed={stats['input_rows']} written={stats['output_rows']}",
                    file=sys.stderr,
                )

    write_summary(summary_path, stats, total_samples, tissue_counts, args)
    print(f"[snaptron-ref] summary: {summary_path}", file=sys.stderr)
    print(f"[snaptron-ref][done] input_rows={stats['input_rows']} output_rows={stats['output_rows']}", file=sys.stderr)
    return stats


def main() -> int:
    csv.field_size_limit(sys.maxsize)
    build_reference(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
