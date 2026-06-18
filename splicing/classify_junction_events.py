#!/usr/bin/env python3
"""Classify spl2 neojunctions with SSNIP-style splice-event labels."""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Tuple


BIN_SIZE = 1_000_000
CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
INPUT_SUFFIX = ".spl2.novel_junctions.tsv"
OUTPUT_SUFFIX = ".spl3.event_annotated.tsv"


class Exon:
    def __init__(self, start: int, end: int, raw_rank: Optional[int]) -> None:
        self.start = start
        self.end = end
        self.raw_rank = raw_rank
        self.rank: int = 0


class Feature:
    def __init__(self, start: int, end: int, lab: str, feature_id: int) -> None:
        self.start = start
        self.end = end
        self.lab = lab
        self.id = feature_id
        prefix = "E" if lab == "exon" else "I"
        self.idx = f"{prefix}{feature_id}"


class Transcript:
    def __init__(
        self,
        transcript_id: str,
        gene_id: str,
        gene_name: str,
        gene_type: str,
        transcript_type: str,
        chrom_key_value: str,
        chrom: str,
        strand: str,
    ) -> None:
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.transcript_type = transcript_type
        self.chrom_key = chrom_key_value
        self.chrom = chrom
        self.strand = strand
        self.exons: List[Exon] = []
        self.features: List[Feature] = []

    @property
    def is_protein_coding(self) -> bool:
        return self.gene_type == "protein_coding" or self.transcript_type == "protein_coding"

    @property
    def start(self) -> int:
        return min(exon.start for exon in self.exons)

    @property
    def end(self) -> int:
        return max(exon.end for exon in self.exons)


class GtfModel:
    def __init__(self, transcripts: Dict[str, Transcript], bins: DefaultDict[Tuple[str, str, int], List[Transcript]]) -> None:
        self.transcripts = transcripts
        self.bins = bins
        self.lookup: Dict[str, List[Transcript]] = defaultdict(list)
        for tx in transcripts.values():
            self.lookup[tx.transcript_id].append(tx)
            self.lookup[strip_version(tx.transcript_id)].append(tx)


class Classification:
    def __init__(
        self,
        event_class: str,
        raw_type: str,
        transcript: Optional[Transcript],
        left_feature: str,
        right_feature: str,
        skipped_exons: Sequence[Feature],
        note: str,
    ) -> None:
        self.event_class = event_class
        self.raw_type = raw_type
        self.transcript = transcript
        self.left_feature = left_feature
        self.right_feature = right_feature
        self.skipped_exons = list(skipped_exons)
        self.note = note


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Classify spl2 neojunction candidates into SSNIP-style classes "
            "such as ES, A3+/A3-, A5+/A5-, junction_in_exon, and junction_in_intron."
        )
    )
    ap.add_argument("--gtf", required=True, help="GENCODE/GTF annotation, optionally .gz")
    ap.add_argument("--root", required=True, help="Root to scan for spl2 *.novel_junctions.tsv files")
    ap.add_argument(
        "--out-dir",
        default="",
        help="Optional root for spl3 outputs. Input paths relative to --root are preserved.",
    )
    ap.add_argument(
        "--sample-filter",
        action="append",
        default=[],
        help="Substring filter applied to spl2 input paths/basenames. Repeatable; any match is kept.",
    )
    ap.add_argument("--input-suffix", default=INPUT_SUFFIX, help="Suffix identifying spl2 input TSVs")
    ap.add_argument("--out-suffix", default=OUTPUT_SUFFIX, help="Suffix for spl3 output TSVs")
    ap.add_argument("--force", action="store_true", help="Overwrite existing spl3 outputs")
    ap.add_argument("--dry-run", action="store_true", help="Report planned files without writing")
    ap.add_argument(
        "--include-noncanonical",
        action="store_true",
        help="Keep non-canonical chromosomes/contigs instead of only 1-22,X,Y.",
    )
    return ap.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8", newline="")


def chrom_key(chrom: str) -> str:
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def ssnip_chrom(chrom: str) -> str:
    key = chrom_key(chrom)
    if key == "MT":
        return "chrM"
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def is_canonical_chrom(chrom: str) -> bool:
    return chrom_key(chrom) in CANONICAL_CHROMS


def strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0] if "." in identifier else identifier


def parse_gtf_attrs(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in raw.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
        elif "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def parse_optional_int(value: str) -> Optional[int]:
    if not value:
        return None
    match = re.search(r"\d+", value)
    if not match:
        return None
    return int(match.group(0))


def finalize_transcript(transcript: Transcript) -> None:
    exons = sorted(transcript.exons, key=lambda exon: (exon.start, exon.end))
    raw_ranks = [exon.raw_rank for exon in exons]
    use_raw = all(rank is not None for rank in raw_ranks) and len(set(raw_ranks)) == len(raw_ranks)

    if use_raw:
        for exon in exons:
            exon.rank = int(exon.raw_rank)  # type: ignore[arg-type]
    else:
        transcript_order = sorted(exons, key=lambda exon: (exon.start, exon.end), reverse=(transcript.strand == "-"))
        for rank, exon in enumerate(transcript_order, start=1):
            exon.rank = rank

    features: List[Feature] = []
    for exon in exons:
        features.append(Feature(exon.start, exon.end, "exon", exon.rank))

    for left, right in zip(exons, exons[1:]):
        intron_start = left.end + 1
        intron_end = right.start - 1
        if intron_start > intron_end:
            continue
        intron_rank = min(left.rank, right.rank)
        features.append(Feature(intron_start, intron_end, "intron", intron_rank))

    transcript.exons = exons
    transcript.features = sorted(features, key=lambda feature: (feature.start, feature.end, feature.lab))


def add_to_bins(bins: DefaultDict[Tuple[str, str, int], List[Transcript]], transcript: Transcript) -> None:
    first_bin = transcript.start // BIN_SIZE
    last_bin = transcript.end // BIN_SIZE
    for bin_no in range(first_bin, last_bin + 1):
        bins[(transcript.chrom_key, transcript.strand, bin_no)].append(transcript)


def build_gtf_model(gtf: Path, include_noncanonical: bool) -> GtfModel:
    transcripts: Dict[str, Transcript] = {}

    with open_text(gtf) as fh:
        for lineno, raw_line in enumerate(fh, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue
            cols = raw_line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue

            chrom = cols[0]
            strand = cols[6]
            if strand not in {"+", "-"}:
                continue
            if not include_noncanonical and not is_canonical_chrom(chrom):
                continue

            attrs = parse_gtf_attrs(cols[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue

            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError as exc:
                raise SystemExit(f"ERROR: non-integer GTF exon coordinate at {gtf}:{lineno}: {exc}") from exc

            tx = transcripts.get(transcript_id)
            if tx is None:
                tx = Transcript(
                    transcript_id=transcript_id,
                    gene_id=attrs.get("gene_id", ""),
                    gene_name=attrs.get("gene_name", attrs.get("gene_id", "")),
                    gene_type=attrs.get("gene_type", attrs.get("gene_biotype", "")),
                    transcript_type=attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
                    chrom_key_value=chrom_key(chrom),
                    chrom=ssnip_chrom(chrom),
                    strand=strand,
                )
                transcripts[transcript_id] = tx

            tx.exons.append(Exon(start, end, parse_optional_int(attrs.get("exon_number", ""))))

    pc_transcripts: Dict[str, Transcript] = {}
    bins: DefaultDict[Tuple[str, str, int], List[Transcript]] = defaultdict(list)
    for transcript_id, tx in transcripts.items():
        if not tx.exons or not tx.is_protein_coding:
            continue
        finalize_transcript(tx)
        pc_transcripts[transcript_id] = tx
        add_to_bins(bins, tx)

    return GtfModel(pc_transcripts, bins)


def discover_inputs(root: Path, filters: Sequence[str], input_suffix: str) -> List[Path]:
    inputs: List[Path] = []
    for path in root.rglob(f"*{input_suffix}"):
        if not path.is_file():
            continue
        if filters and not input_matches_filters(path, filters):
            continue
        inputs.append(path)
    return sorted(inputs)


def input_matches_filters(path: Path, filters: Sequence[str]) -> bool:
    lowered_filters = [filter_value.lower() for filter_value in filters]
    if any(filter_value in str(path).lower() for filter_value in lowered_filters):
        return True

    # spl2 can be written directly under splicing_outdir when STAR root points
    # at one sample folder, so the path may be RNA_TUMOUR.spl2... without Pat101.
    # In that case, use the TSV sample/source columns for patient matching.
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


def output_path_for(input_path: Path, root: Path, out_dir: Optional[Path], input_suffix: str, out_suffix: str) -> Path:
    if input_path.name.endswith(input_suffix):
        out_name = f"{input_path.name[:-len(input_suffix)]}{out_suffix}"
    else:
        out_name = f"{input_path.stem}{out_suffix}"
    if out_dir is None:
        return input_path.with_name(out_name)
    try:
        rel_parent = input_path.parent.relative_to(root)
    except ValueError:
        rel_parent = Path(input_path.parent.name)
    return out_dir / rel_parent / out_name


def parse_transcript_ids(value: str) -> List[str]:
    if not value or value == "NA":
        return []
    return [item for item in value.split(";") if item and item != "NA"]


def parse_junction_coords(row: Dict[str, str]) -> Tuple[str, str, int, int]:
    chrom = row.get("chrom", "")
    strand = row.get("strand", "")
    if row.get("star_intron_start") and row.get("star_intron_end"):
        return chrom, strand, int(row["star_intron_start"]), int(row["star_intron_end"])

    junc_id = row.get("junc_id", "")
    try:
        chrom_part, strand_part, coord_part = junc_id.split(":", 2)
        left_part, right_part = coord_part.split("-", 1)
        return chrom_part, strand_part, int(left_part) + 1, int(right_part)
    except ValueError as exc:
        raise SystemExit(f"ERROR: cannot parse junction coordinates from row junc_id={junc_id!r}") from exc


def query_transcripts(model: GtfModel, row: Dict[str, str], chrom: str, strand: str, intron_start: int, intron_end: int) -> List[Transcript]:
    seen: Dict[str, Transcript] = {}
    for transcript_id in parse_transcript_ids(row.get("transcript_ids", "")):
        for tx in model.lookup.get(transcript_id, []):
            seen[tx.transcript_id] = tx
        for tx in model.lookup.get(strip_version(transcript_id), []):
            seen[tx.transcript_id] = tx

    if not seen:
        ckey = chrom_key(chrom)
        first_bin = intron_start // BIN_SIZE
        last_bin = intron_end // BIN_SIZE
        for bin_no in range(first_bin, last_bin + 1):
            for tx in model.bins.get((ckey, strand, bin_no), []):
                if tx.start <= intron_start and intron_end <= tx.end:
                    seen[tx.transcript_id] = tx

    return sorted(seen.values(), key=lambda tx: (tx.gene_name, tx.transcript_id))


def features_containing(features: Iterable[Feature], pos: int, strict: bool = False) -> List[Feature]:
    if strict:
        return [feature for feature in features if feature.start < pos < feature.end]
    return [feature for feature in features if feature.start <= pos <= feature.end]


def feature_label(feature: Optional[Feature]) -> str:
    if feature is None:
        return "NA"
    return f"{feature.lab}:{feature.idx}:{feature.start}-{feature.end}"


def join_unique(values: Iterable[str]) -> str:
    unique = sorted({value for value in values if value})
    return ";".join(unique) if unique else "NA"


def map_raw_type(raw_type: str) -> str:
    mapping = {
        "A3.gain": "A3+",
        "A3.loss": "A3-",
        "A5.gain": "A5+",
        "A5.loss": "A5-",
        "ES": "ES",
        "JUNC.WITHIN.EXON": "junction_in_exon",
        "JUNC.WITHIN.INTRON": "junction_in_intron",
        "OTHERS": "other",
    }
    return mapping.get(raw_type, "other")


def classify_for_transcript(tx: Transcript, intron_start: int, intron_end: int) -> Classification:
    exons = [feature for feature in tx.features if feature.lab == "exon"]
    introns = [feature for feature in tx.features if feature.lab == "intron"]

    left_boundary_exons = [exon for exon in exons if exon.end + 1 == intron_start]
    right_boundary_exons = [exon for exon in exons if intron_end == exon.start - 1]
    if len(left_boundary_exons) == 1 and len(right_boundary_exons) == 1:
        skipped = [exon for exon in exons if intron_start < exon.start and exon.end < intron_end]
        if skipped:
            left = left_boundary_exons[0]
            right = right_boundary_exons[0]
            note = f"lt.{left.idx}-rt.{right.idx};skipped={','.join(exon.idx for exon in skipped)}"
            return Classification("ES", "ES", tx, feature_label(left), feature_label(right), skipped, note)

    right_start_exons = [exon for exon in exons if intron_end + 1 == exon.start]
    if len(right_start_exons) == 1:
        right = right_start_exons[0]
        left_hits = features_containing(tx.features, intron_start)
        left = left_hits[0] if left_hits else None
        is_gain = (
            left is not None
            and left.lab == "intron"
            and ((tx.strand == "+" and left.id == right.id - 1) or (tx.strand == "-" and left.id == right.id))
        )
        raw_type = ("A5.gain" if is_gain else "A5.loss") if tx.strand == "+" else ("A3.gain" if is_gain else "A3.loss")
        note = f"lt.{left.idx.lower() if left else 'NA'}-rt.{right.idx}"
        return Classification(map_raw_type(raw_type), raw_type, tx, feature_label(left), feature_label(right), [], note)

    left_end_exons = [exon for exon in exons if exon.end == intron_start - 1]
    if len(left_end_exons) == 1:
        left = left_end_exons[0]
        right_hits = features_containing(tx.features, intron_end)
        right = right_hits[0] if right_hits else None
        is_gain = (
            right is not None
            and right.lab == "intron"
            and ((tx.strand == "+" and left.id == right.id) or (tx.strand == "-" and left.id == right.id + 1))
        )
        raw_type = ("A3.gain" if is_gain else "A3.loss") if tx.strand == "+" else ("A5.gain" if is_gain else "A5.loss")
        note = f"lt.{left.idx}-rt.{right.idx.lower() if right else 'NA'}"
        return Classification(map_raw_type(raw_type), raw_type, tx, feature_label(left), feature_label(right), [], note)

    exon_hits = [
        exon for exon in exons if exon.start < intron_start < exon.end and exon.start < intron_end < exon.end
    ]
    if len(exon_hits) == 1:
        exon = exon_hits[0]
        note = f"within.{exon.idx.lower()}"
        return Classification("junction_in_exon", "JUNC.WITHIN.EXON", tx, feature_label(exon), feature_label(exon), [], note)

    intron_hits = [
        intron for intron in introns if intron.start < intron_start < intron.end and intron.start < intron_end < intron.end
    ]
    if len(intron_hits) == 1:
        intron = intron_hits[0]
        note = f"within.{intron.idx.lower()}"
        return Classification(
            "junction_in_intron",
            "JUNC.WITHIN.INTRON",
            tx,
            feature_label(intron),
            feature_label(intron),
            [],
            note,
        )

    left = (features_containing(tx.features, intron_start) or [None])[0]
    right = (features_containing(tx.features, intron_end) or [None])[0]
    note = f"lt.{left.idx.lower() if left else 'NA'}-rt.{right.idx.lower() if right else 'NA'}"
    return Classification("other", "OTHERS", tx, feature_label(left), feature_label(right), [], note)


def classification_priority(classification: Classification) -> Tuple[int, str]:
    priorities = {
        "ES": 0,
        "A3+": 1,
        "A5+": 1,
        "A3-": 2,
        "A5-": 2,
        "junction_in_exon": 3,
        "junction_in_intron": 4,
        "other": 9,
        "no_transcript_match": 10,
    }
    transcript_id = classification.transcript.transcript_id if classification.transcript is not None else ""
    return priorities.get(classification.event_class, 9), transcript_id


def classify_row(model: GtfModel, row: Dict[str, str]) -> Dict[str, str]:
    chrom, strand, intron_start, intron_end = parse_junction_coords(row)
    candidates = query_transcripts(model, row, chrom, strand, intron_start, intron_end)
    classifications = [classify_for_transcript(tx, intron_start, intron_end) for tx in candidates]

    if not classifications:
        best = Classification(
            "no_transcript_match",
            "NO_TRANSCRIPT_MATCH",
            None,
            "NA",
            "NA",
            [],
            "no protein-coding transcript model matched this junction",
        )
        classifications = [best]
    else:
        best = sorted(classifications, key=classification_priority)[0]

    matched = [classification.transcript for classification in classifications if classification.transcript is not None]
    skipped_ids = [exon.idx for exon in best.skipped_exons]
    all_pairs = []
    for classification in classifications:
        transcript_id = classification.transcript.transcript_id if classification.transcript is not None else "NA"
        all_pairs.append(f"{transcript_id}:{classification.event_class}")

    return {
        "ssnip_event_class": best.event_class,
        "ssnip_event_raw": best.raw_type,
        "matched_transcript_id": best.transcript.transcript_id if best.transcript is not None else "NA",
        "matched_gene_id": best.transcript.gene_id if best.transcript is not None else "NA",
        "matched_gene_name": best.transcript.gene_name if best.transcript is not None else "NA",
        "matched_transcript_count": str(len(matched)),
        "left_feature": best.left_feature,
        "right_feature": best.right_feature,
        "skipped_exon_count": str(len(best.skipped_exons)),
        "skipped_exon_ids": ";".join(skipped_ids) if skipped_ids else "NA",
        "all_ssnip_event_classes": join_unique(classification.event_class for classification in classifications),
        "all_matched_transcripts": ";".join(all_pairs) if all_pairs else "NA",
        "classification_note": best.note,
    }


def process_file(model: GtfModel, input_path: Path, out_path: Path) -> Dict[str, int]:
    stats: DefaultDict[str, int] = defaultdict(int)
    with input_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: missing TSV header in {input_path}")
        added_columns = [
            "ssnip_event_class",
            "ssnip_event_raw",
            "matched_transcript_id",
            "matched_gene_id",
            "matched_gene_name",
            "matched_transcript_count",
            "left_feature",
            "right_feature",
            "skipped_exon_count",
            "skipped_exon_ids",
            "all_ssnip_event_classes",
            "all_matched_transcripts",
            "classification_note",
        ]
        fieldnames = list(reader.fieldnames) + [col for col in added_columns if col not in reader.fieldnames]

        rows: List[Dict[str, str]] = []
        for row in reader:
            stats["input_rows"] += 1
            annotations = classify_row(model, row)
            stats[f"class_{annotations['ssnip_event_class']}"] += 1
            row.update(annotations)
            rows.append(row)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    return dict(stats)


def print_stats(path: Path, stats: Dict[str, int]) -> None:
    class_fields = sorted(field for field in stats if field.startswith("class_"))
    fields = ["input_rows"] + class_fields
    details = " ".join(f"{field}={stats.get(field, 0)}" for field in fields)
    print(f"[spl3] {path}: {details}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    gtf = Path(args.gtf)
    root = Path(args.root)
    out_dir = Path(args.out_dir) if args.out_dir else None
    if not gtf.exists():
        raise SystemExit(f"ERROR: GTF does not exist: {gtf}")
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"ERROR: spl2 root does not exist or is not a directory: {root}")

    print(f"[spl3] loading GTF: {gtf}", file=sys.stderr)
    model = build_gtf_model(gtf, include_noncanonical=args.include_noncanonical)
    print(f"[spl3] GTF protein-coding transcripts={len(model.transcripts)}", file=sys.stderr)

    inputs = discover_inputs(root, args.sample_filter, args.input_suffix)
    if not inputs:
        print(f"[spl3] no spl2 TSV files found under {root}", file=sys.stderr)
        return 0

    written = 0
    skipped = 0
    for input_path in inputs:
        out_path = output_path_for(input_path, root, out_dir, args.input_suffix, args.out_suffix)
        if out_path.exists() and not args.force:
            print(f"[spl3][skip] output exists: {out_path}", file=sys.stderr)
            skipped += 1
            continue
        print(f"[spl3] {out_path} <- {input_path}", file=sys.stderr)
        if args.dry_run:
            written += 1
            continue
        stats = process_file(model, input_path, out_path)
        print_stats(out_path, stats)
        written += 1

    print(f"[spl3][done] written={written} skipped={skipped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
