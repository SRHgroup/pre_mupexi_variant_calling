#!/usr/bin/env python3
"""Filter spl3 neojunctions against a normal splice-junction reference."""

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Optional, Sequence, Tuple


INPUT_SUFFIX = ".spl3.event_annotated.tsv"
OUTPUT_SUFFIX = ".spl3.5.cancer_unique.tsv"
FILTERED_SUFFIX = ".spl3.5.normal_present.tsv"
SUMMARY_SUFFIX = ".spl3.5.normal_filter_summary.tsv"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Filter spl3 SSNIP-style neojunctions against a compact normal "
            "splice-junction reference before spl4 sequence reconstruction."
        )
    )
    ap.add_argument("--root", required=True, help="Root to scan for spl3 TSV files")
    ap.add_argument("--normal-ref", required=True, help="Normal junction reference TSV/TSV.GZ")
    ap.add_argument(
        "--out-dir",
        default="",
        help="Optional root for spl3.5 outputs. Input paths relative to --root are preserved.",
    )
    ap.add_argument(
        "--output-subdir",
        default="",
        help="Optional first-level folder under --out-dir for patient-scoped outputs.",
    )
    ap.add_argument("--sample-filter", action="append", default=[], help="Substring filter for input paths")
    ap.add_argument("--input-suffix", default=INPUT_SUFFIX, help="Suffix identifying spl3 input TSVs")
    ap.add_argument("--out-suffix", default=OUTPUT_SUFFIX, help="Suffix for cancer-unique output TSVs")
    ap.add_argument("--filtered-suffix", default=FILTERED_SUFFIX, help="Suffix for normal-present audit TSVs")
    ap.add_argument("--summary-suffix", default=SUMMARY_SUFFIX, help="Suffix for per-file summary TSVs")
    ap.add_argument(
        "--max-normal-prevalence",
        type=float,
        default=0.01,
        help="Keep junctions only when normal prevalence is below this value. Default: 0.01.",
    )
    ap.add_argument(
        "--normal-total-samples",
        type=float,
        default=0.0,
        help="Total normal samples used to calculate prevalence when the reference has counts but no prevalence column.",
    )
    ap.add_argument(
        "--max-normal-sample-count",
        type=float,
        default=None,
        help="Optional count threshold. If prevalence is unavailable, filter rows above this count.",
    )
    ap.add_argument("--ignore-strand", action="store_true", help="Match normal reference junctions without strand")
    ap.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    ap.add_argument("--dry-run", action="store_true", help="Report planned files without writing")
    return ap.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r", encoding="utf-8", newline="")


def chrom_key(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def parse_float(value: str) -> Optional[float]:
    if value is None:
        return None
    value = str(value).strip()
    if not value or value.upper() in {"NA", "NAN", "."}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int_string(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.6g}"


def first_present(row: Dict[str, str], columns: Sequence[str]) -> str:
    lower = {key.lower(): key for key in row}
    for column in columns:
        key = lower.get(column.lower())
        if key is not None and row.get(key, "") != "":
            return row[key]
    return ""


def parse_boundary(value: str) -> Optional[int]:
    parsed = parse_float(value)
    if parsed is None:
        return None
    return int(parsed)


def row_key(row: Dict[str, str], ignore_strand: bool) -> Optional[Tuple[str, int, int, str]]:
    chrom = first_present(row, ["chrom", "chromosome", "chr"])
    if not chrom:
        return None

    left = first_present(row, ["left_boundary", "start", "junction_start", "intron_start", "star_intron_start"])
    right = first_present(row, ["right_boundary", "end", "junction_end", "intron_end", "star_intron_end"])
    left_i = parse_boundary(left)
    right_i = parse_boundary(right)
    if left_i is None or right_i is None:
        return None

    if left_i > right_i:
        left_i, right_i = right_i, left_i

    strand = "." if ignore_strand else first_present(row, ["strand"])
    if not strand or strand in {"?", "NA"}:
        strand = "."
    return chrom_key(chrom), left_i, right_i, strand


class NormalHit:
    def __init__(self) -> None:
        self.sample_count: Optional[float] = None
        self.total_reads: Optional[float] = None
        self.prevalence: Optional[float] = None
        self.sources: List[str] = []
        self.tissue_sample_counts: Counter = Counter()
        self.tissue_read_counts: Counter = Counter()
        self.broad_tissue_sample_counts: Counter = Counter()
        self.broad_tissue_read_counts: Counter = Counter()

    def add(self, row: Dict[str, str], total_samples: float) -> None:
        sample_count = parse_float(
            first_present(row, ["normal_sample_count", "gtex_sample_count", "sample_count", "samples_count"])
        )
        total_reads = parse_float(
            first_present(row, ["normal_total_reads", "gtex_total_reads", "coverage_sum", "total_reads", "read_count"])
        )
        prevalence = parse_float(
            first_present(row, ["normal_prevalence", "gtex_prevalence", "prevalence", "psr", "psr_gtex"])
        )
        if prevalence is None and sample_count is not None and total_samples > 0:
            prevalence = sample_count / total_samples

        self.sample_count = max_optional(self.sample_count, sample_count)
        self.total_reads = sum_optional(self.total_reads, total_reads)
        self.prevalence = max_optional(self.prevalence, prevalence)
        source = first_present(row, ["source", "normal_source", "gtex_source", "dataset"])
        if source and source not in self.sources:
            self.sources.append(source)
        self.tissue_sample_counts.update(parse_named_counts(first_present(row, ["normal_tissue_sample_counts"])))
        self.tissue_read_counts.update(parse_named_counts(first_present(row, ["normal_tissue_read_counts"])))
        self.broad_tissue_sample_counts.update(parse_named_counts(first_present(row, ["normal_broad_tissue_sample_counts"])))
        self.broad_tissue_read_counts.update(parse_named_counts(first_present(row, ["normal_broad_tissue_read_counts"])))


def parse_named_counts(value: str) -> Counter:
    counts: Counter = Counter()
    if not value or value.upper() in {"NA", "NAN", "."}:
        return counts
    for item in value.split(";"):
        if not item or "=" not in item:
            continue
        label, raw_count = item.rsplit("=", 1)
        parsed = parse_float(raw_count)
        if parsed is not None:
            counts[label] += parsed
    return counts


def format_named_counts(counter: Counter) -> str:
    if not counter:
        return "NA"
    parts = []
    for label, value in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
        if float(value).is_integer():
            value_s = str(int(value))
        else:
            value_s = f"{value:.6g}"
        parts.append(f"{label}={value_s}")
    return ";".join(parts)


def label_count(counter: Counter) -> str:
    return str(len(counter))


def max_optional(left: Optional[float], right: Optional[float]) -> Optional[float]:
    if left is None:
        return right
    if right is None:
        return left
    return max(left, right)


def sum_optional(left: Optional[float], right: Optional[float]) -> Optional[float]:
    if left is None:
        return right
    if right is None:
        return left
    return left + right


def load_normal_reference(path: Path, total_samples: float, ignore_strand: bool) -> Dict[Tuple[str, int, int, str], NormalHit]:
    hits: Dict[Tuple[str, int, int, str], NormalHit] = {}
    skipped = 0
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: normal reference has no header: {path}")
        for row in reader:
            key = row_key(row, ignore_strand)
            if key is None:
                skipped += 1
                continue
            hit = hits.setdefault(key, NormalHit())
            hit.add(row, total_samples)
    print(f"[spl3.5] normal reference rows_loaded={len(hits)} skipped={skipped}", file=sys.stderr)
    return hits


def input_matches_filters(path: Path, filters: Sequence[str]) -> bool:
    if not filters:
        return True
    lowered_filters = [filter_value.lower() for filter_value in filters]
    if any(filter_value in str(path).lower() for filter_value in lowered_filters):
        return True
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                haystack = " ".join(
                    row.get(column, "")
                    for column in ("sample", "source_sj", "junc_id", "gene_names", "transcript_ids")
                ).lower()
                if any(filter_value in haystack for filter_value in lowered_filters):
                    return True
    except OSError:
        return False
    return False


def infer_input_sample_label(path: Path, input_suffix: str) -> str:
    name = path.name
    if name.endswith(input_suffix):
        name = name[: -len(input_suffix)]
    return name


def input_source_priority(path: Path, input_suffix: str) -> Tuple[int, int, str]:
    sample = infer_input_sample_label(path, input_suffix)
    parent = path.parent.name
    score = 1 if parent == sample else 0
    return score, len(str(path)), str(path)


def discover_inputs(root: Path, filters: Sequence[str], input_suffix: str) -> Tuple[List[Path], Dict[str, List[Path]]]:
    grouped: DefaultDict[str, List[Path]] = defaultdict(list)
    for path in root.rglob(f"*{input_suffix}"):
        if path.is_file() and input_matches_filters(path, filters):
            grouped[infer_input_sample_label(path, input_suffix)].append(path)

    selected: List[Path] = []
    duplicates: Dict[str, List[Path]] = {}
    for sample, candidates in grouped.items():
        ranked = sorted(candidates, key=lambda candidate: input_source_priority(candidate, input_suffix), reverse=True)
        selected.append(ranked[0])
        if len(ranked) > 1:
            duplicates[sample] = ranked
    return sorted(selected), duplicates


def scoped_rel_parent(input_parent: Path, root: Path, output_subdir: str) -> Path:
    try:
        rel_parent = input_parent.relative_to(root)
    except ValueError:
        rel_parent = Path(input_parent.name)
    if output_subdir:
        parts = rel_parent.parts
        if not parts or parts[0] != output_subdir:
            return Path(output_subdir)
    return rel_parent


def output_path_for(
    input_path: Path,
    root: Path,
    out_dir: Optional[Path],
    input_suffix: str,
    out_suffix: str,
    output_subdir: str,
) -> Path:
    if input_path.name.endswith(input_suffix):
        out_name = f"{input_path.name[:-len(input_suffix)]}{out_suffix}"
    else:
        out_name = f"{input_path.stem}{out_suffix}"
    if out_dir is None:
        return input_path.with_name(out_name)
    rel_parent = scoped_rel_parent(input_path.parent, root, output_subdir)
    return out_dir / rel_parent / out_name


def is_normal_positive(hit: NormalHit, args: argparse.Namespace) -> Tuple[bool, str]:
    if hit.prevalence is not None:
        if hit.prevalence >= args.max_normal_prevalence:
            return True, f"normal_prevalence>={args.max_normal_prevalence:g}"
        return False, f"normal_prevalence<{args.max_normal_prevalence:g}"
    if hit.sample_count is not None and args.normal_total_samples > 0:
        prevalence = hit.sample_count / args.normal_total_samples
        if prevalence >= args.max_normal_prevalence:
            return True, f"normal_count/{args.normal_total_samples:g}>={args.max_normal_prevalence:g}"
        return False, f"normal_count/{args.normal_total_samples:g}<{args.max_normal_prevalence:g}"
    if hit.sample_count is not None and args.max_normal_sample_count is not None:
        if hit.sample_count > args.max_normal_sample_count:
            return True, f"normal_sample_count>{args.max_normal_sample_count:g}"
        return False, f"normal_sample_count<={args.max_normal_sample_count:g}"
    if hit.sample_count is not None:
        if hit.sample_count > 0:
            return True, "normal_sample_count>0"
        return False, "normal_sample_count=0"
    return True, "normal_reference_match"


def annotate_row(row: Dict[str, str], hit: Optional[NormalHit], args: argparse.Namespace) -> Tuple[Dict[str, str], bool]:
    row = dict(row)
    if hit is None:
        row.update(
            {
                "normal_ref_match": "0",
                "normal_status": "absent",
                "normal_filter_reason": "not_in_normal_reference",
                "normal_sample_count": "0",
                "normal_total_reads": "0",
                "normal_prevalence": "0",
                "normal_source": "NA",
                "normal_tissue_count": "0",
                "normal_tissue_sample_counts": "NA",
                "normal_tissue_read_counts": "NA",
                "normal_broad_tissue_count": "0",
                "normal_broad_tissue_sample_counts": "NA",
                "normal_broad_tissue_read_counts": "NA",
            }
        )
        return row, True

    positive, reason = is_normal_positive(hit, args)
    row.update(
        {
            "normal_ref_match": "1",
            "normal_status": "filtered_normal_present" if positive else "present_below_threshold",
            "normal_filter_reason": reason,
            "normal_sample_count": parse_int_string(hit.sample_count),
            "normal_total_reads": parse_int_string(hit.total_reads),
            "normal_prevalence": "NA" if hit.prevalence is None else f"{hit.prevalence:.8g}",
            "normal_source": ";".join(hit.sources) if hit.sources else "normal_reference",
            "normal_tissue_count": label_count(hit.tissue_sample_counts),
            "normal_tissue_sample_counts": format_named_counts(hit.tissue_sample_counts),
            "normal_tissue_read_counts": format_named_counts(hit.tissue_read_counts),
            "normal_broad_tissue_count": label_count(hit.broad_tissue_sample_counts),
            "normal_broad_tissue_sample_counts": format_named_counts(hit.broad_tissue_sample_counts),
            "normal_broad_tissue_read_counts": format_named_counts(hit.broad_tissue_read_counts),
        }
    )
    return row, not positive


def append_fields(fieldnames: Sequence[str], extra: Sequence[str]) -> List[str]:
    fields = list(fieldnames)
    for field in extra:
        if field not in fields:
            fields.append(field)
    return fields


def write_rows(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def update_tissue_summary_stats(stats, annotated: Dict[str, str], keep: bool) -> None:
    if annotated.get("normal_ref_match") != "1":
        return
    for label in parse_named_counts(annotated.get("normal_tissue_sample_counts", "")).keys():
        stats[f"normal_tissue_match:{label}"] += 1
        if not keep:
            stats[f"normal_tissue_filtered:{label}"] += 1
    for label in parse_named_counts(annotated.get("normal_broad_tissue_sample_counts", "")).keys():
        stats[f"normal_broad_tissue_match:{label}"] += 1
        if not keep:
            stats[f"normal_broad_tissue_filtered:{label}"] += 1


def process_file(
    input_path: Path,
    output_path: Path,
    filtered_path: Path,
    summary_path: Path,
    normal_hits: Dict[Tuple[str, int, int, str], NormalHit],
    args: argparse.Namespace,
) -> Dict[str, int]:
    kept_rows: List[Dict[str, str]] = []
    filtered_rows: List[Dict[str, str]] = []
    stats = defaultdict(int)
    with input_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: input has no header: {input_path}")
        fields = append_fields(
            reader.fieldnames,
            [
                "normal_ref_match",
                "normal_status",
                "normal_filter_reason",
                "normal_sample_count",
                "normal_total_reads",
                "normal_prevalence",
                "normal_source",
                "normal_tissue_count",
                "normal_tissue_sample_counts",
                "normal_tissue_read_counts",
                "normal_broad_tissue_count",
                "normal_broad_tissue_sample_counts",
                "normal_broad_tissue_read_counts",
            ],
        )
        for row in reader:
            stats["input_rows"] += 1
            key = row_key(row, args.ignore_strand)
            hit = normal_hits.get(key) if key is not None else None
            annotated, keep = annotate_row(row, hit, args)
            if keep:
                kept_rows.append(annotated)
                stats["kept_rows"] += 1
            else:
                filtered_rows.append(annotated)
                stats["filtered_rows"] += 1
            stats[f"normal_status_{annotated['normal_status']}"] += 1
            update_tissue_summary_stats(stats, annotated, keep)

    write_rows(output_path, fields, kept_rows)
    write_rows(filtered_path, fields, filtered_rows)
    summary_fields = ["input", "output", "normal_present_output", "metric", "value"]
    summary_rows = [
        {
            "input": str(input_path),
            "output": str(output_path),
            "normal_present_output": str(filtered_path),
            "metric": metric,
            "value": str(value),
        }
        for metric, value in sorted(stats.items())
    ]
    write_rows(summary_path, summary_fields, summary_rows)
    return dict(stats)


def print_stats(path: Path, stats: Dict[str, int]) -> None:
    fields = ["input_rows", "kept_rows", "filtered_rows"]
    details = " ".join(f"{field}={stats.get(field, 0)}" for field in fields)
    print(f"[spl3.5] {path}: {details}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    csv.field_size_limit(sys.maxsize)
    root = Path(args.root)
    normal_ref = Path(args.normal_ref)
    out_dir = Path(args.out_dir) if args.out_dir else None
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"ERROR: spl3 root does not exist or is not a directory: {root}")
    if not normal_ref.exists():
        raise SystemExit(f"ERROR: normal reference does not exist: {normal_ref}")

    print(f"[spl3.5] loading normal reference: {normal_ref}", file=sys.stderr)
    normal_hits = load_normal_reference(normal_ref, args.normal_total_samples, args.ignore_strand)

    inputs, duplicate_inputs = discover_inputs(root, args.sample_filter, args.input_suffix)
    if not inputs:
        print(f"[spl3.5] no spl3 TSV files found under {root}", file=sys.stderr)
        return 0
    for sample, ranked in sorted(duplicate_inputs.items()):
        skipped = ", ".join(str(path) for path in ranked[1:])
        print(f"[spl3.5][dedupe] sample={sample} keeping {ranked[0]} skipped={skipped}", file=sys.stderr)

    written = 0
    skipped = 0
    for input_path in inputs:
        out_path = output_path_for(input_path, root, out_dir, args.input_suffix, args.out_suffix, args.output_subdir)
        filtered_path = output_path_for(input_path, root, out_dir, args.input_suffix, args.filtered_suffix, args.output_subdir)
        summary_path = output_path_for(input_path, root, out_dir, args.input_suffix, args.summary_suffix, args.output_subdir)
        outputs_exist = out_path.exists() and filtered_path.exists() and summary_path.exists()
        if outputs_exist and not args.force:
            print(f"[spl3.5][skip] outputs exist: {out_path}", file=sys.stderr)
            skipped += 1
            continue
        print(f"[spl3.5] {out_path} <- {input_path}", file=sys.stderr)
        if args.dry_run:
            written += 1
            continue
        stats = process_file(input_path, out_path, filtered_path, summary_path, normal_hits, args)
        print_stats(out_path, stats)
        written += 1

    print(f"[spl3.5][done] written={written} skipped={skipped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
