#!/usr/bin/env python3
"""Call SSNIP-style tumor-supported novel splice junction candidates."""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Set, Tuple


SHARD_RE = re.compile(r"^.+\.\d{4}\.SJ\.out\.tab(?:\.gz)?$")
SJ_RE = re.compile(r"^.+\.SJ\.out\.tab(?:\.gz)?$")
BIN_SIZE = 1_000_000

CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}
STAR_MOTIFS = {
    0: "noncanonical",
    1: "GT/AG",
    2: "CT/AC",
    3: "GC/AG",
    4: "CT/GC",
    5: "AT/AC",
    6: "GT/AT",
}


class TranscriptInfo:
    def __init__(
        self,
        transcript_id: str,
        gene_id: str,
        gene_name: str,
        gene_type: str,
        transcript_type: str,
        chrom_key: str,
        chrom: str,
        strand: str,
    ) -> None:
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.transcript_type = transcript_type
        self.chrom_key = chrom_key
        self.chrom = chrom
        self.strand = strand
        self.exons: List[Tuple[int, int]] = []

    @property
    def start(self) -> int:
        return min(start for start, _ in self.exons)

    @property
    def end(self) -> int:
        return max(end for _, end in self.exons)

    @property
    def is_protein_coding(self) -> bool:
        return self.gene_type == "protein_coding" or self.transcript_type == "protein_coding"


class GtfIndex:
    def __init__(
        self,
        known_junctions: Set[Tuple[str, str, int, int]],
        left_edges: DefaultDict[Tuple[str, str], Set[int]],
        right_edges: DefaultDict[Tuple[str, str], Set[int]],
        pc_bins: DefaultDict[Tuple[str, str, int], List[TranscriptInfo]],
        transcript_count: int,
        protein_coding_transcript_count: int,
        known_junction_count: int,
    ) -> None:
        self.known_junctions = known_junctions
        self.left_edges = left_edges
        self.right_edges = right_edges
        self.pc_bins = pc_bins
        self.transcript_count = transcript_count
        self.protein_coding_transcript_count = protein_coding_transcript_count
        self.known_junction_count = known_junction_count


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Reimplement the early SSNIP neojunction-calling logic: build a "
            "GTF-derived known splice-junction set, then report STAR junctions "
            "that are tumor-supported, non-annotated, and overlapping "
            "protein-coding transcript spans."
        )
    )
    ap.add_argument("--gtf", required=True, help="GENCODE/GTF annotation, optionally .gz")
    ap.add_argument("--star-root", required=True, help="Root to scan for merged STAR SJ.out.tab files")
    ap.add_argument(
        "--sample-filter",
        action="append",
        default=[],
        help="Substring filter applied to merged SJ paths/basenames. Repeatable; any match is kept.",
    )
    ap.add_argument(
        "--min-unique-reads",
        type=int,
        default=10,
        help="Minimum STAR uniquely mapped junction reads. SSNIP used 10.",
    )
    ap.add_argument(
        "--out-suffix",
        default=".spl2.novel_junctions.tsv",
        help="Suffix added to merged SJ basename for the output TSV.",
    )
    ap.add_argument(
        "--out-dir",
        default="",
        help=(
            "Optional root directory for spl2 TSV outputs. When set, the "
            "input path relative to --star-root is preserved under this root."
        ),
    )
    ap.add_argument(
        "--output-subdir",
        default="",
        help=(
            "Optional first-level folder under --out-dir for patient-scoped "
            "outputs, e.g. Pat101_RNA_TUMOR. If the input is already under "
            "that folder relative to the scan root, it is not added twice."
        ),
    )
    ap.add_argument("--force", action="store_true", help="Overwrite existing spl2 outputs")
    ap.add_argument("--dry-run", action="store_true", help="Report planned files without writing outputs")
    ap.add_argument(
        "--include-noncanonical",
        action="store_true",
        help="Keep non-canonical chromosomes/contigs instead of only 1-22,X,Y.",
    )
    ap.add_argument(
        "--keep-non-protein-coding",
        action="store_true",
        help="Keep novel junctions even if they do not overlap a protein-coding transcript span.",
    )
    return ap.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8")


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


def chrom_sort_key(chrom: str) -> Tuple[int, object]:
    key = chrom_key(chrom)
    if key.isdigit():
        return (0, int(key))
    if key == "X":
        return (0, 23)
    if key == "Y":
        return (0, 24)
    if key == "MT":
        return (0, 25)
    return (1, chrom)


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


def add_to_bins(
    bins: DefaultDict[Tuple[str, str, int], List[TranscriptInfo]],
    transcript: TranscriptInfo,
) -> None:
    first_bin = transcript.start // BIN_SIZE
    last_bin = transcript.end // BIN_SIZE
    for bin_no in range(first_bin, last_bin + 1):
        bins[(transcript.chrom_key, transcript.strand, bin_no)].append(transcript)


def build_gtf_index(gtf: Path, include_noncanonical: bool) -> GtfIndex:
    transcripts: Dict[str, TranscriptInfo] = {}

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

            info = transcripts.get(transcript_id)
            if info is None:
                info = TranscriptInfo(
                    transcript_id=transcript_id,
                    gene_id=attrs.get("gene_id", ""),
                    gene_name=attrs.get("gene_name", attrs.get("gene_id", "")),
                    gene_type=attrs.get("gene_type", attrs.get("gene_biotype", "")),
                    transcript_type=attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
                    chrom_key=chrom_key(chrom),
                    chrom=ssnip_chrom(chrom),
                    strand=strand,
                )
                transcripts[transcript_id] = info

            info.exons.append((start, end))

    known_junctions: Set[Tuple[str, str, int, int]] = set()
    left_edges: DefaultDict[Tuple[str, str], Set[int]] = defaultdict(set)
    right_edges: DefaultDict[Tuple[str, str], Set[int]] = defaultdict(set)
    pc_bins: DefaultDict[Tuple[str, str, int], List[TranscriptInfo]] = defaultdict(list)
    pc_count = 0

    for transcript in transcripts.values():
        if len(transcript.exons) >= 2:
            exons = sorted(transcript.exons)
            for left_exon, right_exon in zip(exons, exons[1:]):
                intron_start = left_exon[1] + 1
                intron_end = right_exon[0] - 1
                if intron_start > intron_end:
                    continue
                left_boundary = intron_start - 1
                right_boundary = intron_end
                key = (transcript.chrom_key, transcript.strand, left_boundary, right_boundary)
                known_junctions.add(key)
                left_edges[(transcript.chrom_key, transcript.strand)].add(left_boundary)
                right_edges[(transcript.chrom_key, transcript.strand)].add(right_boundary)

        if transcript.is_protein_coding and transcript.exons:
            pc_count += 1
            add_to_bins(pc_bins, transcript)

    return GtfIndex(
        known_junctions=known_junctions,
        left_edges=left_edges,
        right_edges=right_edges,
        pc_bins=pc_bins,
        transcript_count=len(transcripts),
        protein_coding_transcript_count=pc_count,
        known_junction_count=len(known_junctions),
    )


def sj_source_priority(path: Path) -> Tuple[int, int, str]:
    sample = infer_sample_label(path).lower()
    base = sj_base(path).lower()
    parent = path.parent.name.lower()
    score = 0
    if base == sample:
        score += 100
    if parent == sample:
        score += 50
    if base.startswith(sample):
        score += 20
    if parent.startswith(sample):
        score += 10
    try:
        size = path.stat().st_size
    except OSError:
        size = 0
    return (score, size, str(path))


def discover_sj_files(root: Path, filters: Sequence[str]) -> Tuple[List[Path], Dict[str, List[Path]]]:
    grouped: DefaultDict[str, List[Path]] = defaultdict(list)
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if SHARD_RE.match(path.name):
            continue
        if not SJ_RE.match(path.name):
            continue
        if filters and not any(f.lower() in str(path).lower() for f in filters):
            continue
        grouped[infer_sample_label(path)].append(path)

    selected: List[Path] = []
    duplicates: Dict[str, List[Path]] = {}
    for sample, candidates in grouped.items():
        ranked = sorted(candidates, key=sj_source_priority, reverse=True)
        selected.append(ranked[0])
        if len(ranked) > 1:
            duplicates[sample] = ranked
    return sorted(selected), duplicates


def sj_base(path: Path) -> str:
    name = path.name
    for suffix in (".SJ.out.tab.gz", ".SJ.out.tab"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


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
    sj_path: Path,
    star_root: Path,
    out_dir: Optional[Path],
    out_suffix: str,
    output_subdir: str,
) -> Path:
    out_name = f"{sj_base(sj_path)}{out_suffix}"
    if not out_dir:
        return sj_path.with_name(out_name)
    rel_parent = scoped_rel_parent(sj_path.parent, star_root, output_subdir)
    return out_dir / rel_parent / out_name


def infer_sample_label(sj_path: Path) -> str:
    base = sj_base(sj_path)
    parent = sj_path.parent.name
    grandparent = sj_path.parent.parent.name if sj_path.parent.parent != sj_path.parent else ""
    if (
        grandparent
        and parent == base
        and not base.startswith(f"{grandparent}_")
        and grandparent.lower() not in {"star", "reports"}
    ):
        return f"{grandparent}_{base}"
    return base


def star_strand(value: str) -> str:
    if value == "1":
        return "+"
    if value == "2":
        return "-"
    return "undefined"


def query_pc_transcripts(
    index: GtfIndex,
    chrom: str,
    strand: str,
    intron_start: int,
    intron_end: int,
) -> List[TranscriptInfo]:
    first_bin = intron_start // BIN_SIZE
    last_bin = intron_end // BIN_SIZE
    seen: Set[str] = set()
    hits: List[TranscriptInfo] = []
    for bin_no in range(first_bin, last_bin + 1):
        for transcript in index.pc_bins.get((chrom, strand, bin_no), []):
            if transcript.transcript_id in seen:
                continue
            seen.add(transcript.transcript_id)
            if transcript.start < intron_start and intron_end < transcript.end:
                hits.append(transcript)
    return sorted(hits, key=lambda tx: (tx.gene_name, tx.transcript_id))


def join_unique(values: Iterable[str]) -> str:
    unique = sorted({value for value in values if value})
    return ";".join(unique) if unique else "NA"


def event_class(strand: str, left_known: bool, right_known: bool) -> str:
    if strand == "+":
        donor_known = left_known
        acceptor_known = right_known
    else:
        donor_known = right_known
        acceptor_known = left_known

    if donor_known and acceptor_known:
        return "novel_pair_known_edges"
    if donor_known:
        return "alternative_3prime_acceptor"
    if acceptor_known:
        return "alternative_5prime_donor"
    return "novel_donor_and_acceptor"


def output_header() -> List[str]:
    return [
        "sample",
        "source_sj",
        "junc_id",
        "chrom",
        "strand",
        "star_intron_start",
        "star_intron_end",
        "left_boundary",
        "right_boundary",
        "donor_boundary",
        "acceptor_boundary",
        "motif_code",
        "motif",
        "star_annotated",
        "gtf_annotated",
        "unique_reads",
        "multimap_reads",
        "total_reads",
        "max_splice_overhang",
        "left_edge_known",
        "right_edge_known",
        "event_class",
        "gene_ids",
        "gene_names",
        "transcript_ids",
        "protein_coding_overlap_count",
        "normal_status",
        "normal_sample_count",
        "normal_total_reads",
        "ssnip_count_pass",
    ]


def process_sj_file(
    sj_path: Path,
    out_path: Path,
    index: GtfIndex,
    min_unique_reads: int,
    include_noncanonical: bool,
    keep_non_protein_coding: bool,
) -> Dict[str, int]:
    sample = infer_sample_label(sj_path)
    stats = defaultdict(int)
    rows: List[List[str]] = []

    with open_text(sj_path) as fh:
        for lineno, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            stats["input_junctions"] += 1
            cols = line.split("\t")
            if len(cols) != 9:
                raise SystemExit(f"ERROR: expected 9 STAR columns in {sj_path}:{lineno}, got {len(cols)}")

            chrom_raw = cols[0]
            if not include_noncanonical and not is_canonical_chrom(chrom_raw):
                stats["skip_noncanonical_chrom"] += 1
                continue

            strand = star_strand(cols[3])
            if strand == "undefined":
                stats["skip_undefined_strand"] += 1
                continue

            try:
                intron_start = int(cols[1])
                intron_end = int(cols[2])
                motif_code = int(cols[4])
                star_annotated = int(cols[5])
                unique_reads = int(cols[6])
                multimap_reads = int(cols[7])
                max_overhang = int(cols[8])
            except ValueError as exc:
                raise SystemExit(f"ERROR: non-integer STAR field in {sj_path}:{lineno}: {exc}") from exc

            if unique_reads < min_unique_reads:
                stats["skip_low_unique_reads"] += 1
                continue

            chrom = ssnip_chrom(chrom_raw)
            ckey = chrom_key(chrom_raw)
            left_boundary = intron_start - 1
            right_boundary = intron_end
            key = (ckey, strand, left_boundary, right_boundary)
            gtf_annotated = key in index.known_junctions
            if gtf_annotated or star_annotated == 1:
                stats["skip_annotated"] += 1
                continue

            pc_hits = query_pc_transcripts(index, ckey, strand, intron_start, intron_end)
            if not pc_hits and not keep_non_protein_coding:
                stats["skip_no_protein_coding_overlap"] += 1
                continue

            left_known = left_boundary in index.left_edges.get((ckey, strand), set())
            right_known = right_boundary in index.right_edges.get((ckey, strand), set())
            donor_boundary = left_boundary if strand == "+" else right_boundary
            acceptor_boundary = right_boundary if strand == "+" else left_boundary
            junc_id = f"{chrom}:{strand}:{left_boundary}-{right_boundary}"
            total_reads = unique_reads + multimap_reads

            rows.append(
                [
                    sample,
                    str(sj_path),
                    junc_id,
                    chrom,
                    strand,
                    str(intron_start),
                    str(intron_end),
                    str(left_boundary),
                    str(right_boundary),
                    str(donor_boundary),
                    str(acceptor_boundary),
                    str(motif_code),
                    STAR_MOTIFS.get(motif_code, "unknown"),
                    str(star_annotated),
                    "1" if gtf_annotated else "0",
                    str(unique_reads),
                    str(multimap_reads),
                    str(total_reads),
                    str(max_overhang),
                    "1" if left_known else "0",
                    "1" if right_known else "0",
                    event_class(strand, left_known, right_known),
                    join_unique(tx.gene_id for tx in pc_hits),
                    join_unique(tx.gene_name for tx in pc_hits),
                    join_unique(tx.transcript_id for tx in pc_hits),
                    str(len(pc_hits)),
                    "not_checked",
                    "NA",
                    "NA",
                    "1",
                ]
            )
            stats["novel_candidates"] += 1

    rows.sort(
        key=lambda row: (
            chrom_sort_key(row[3]),
            int(row[5]),
            int(row[6]),
            row[4],
            -int(row[15]),
        )
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write("\t".join(output_header()))
        out.write("\n")
        for row in rows:
            out.write("\t".join(row))
            out.write("\n")

    return dict(stats)


def print_stats(path: Path, stats: Dict[str, int]) -> None:
    fields = [
        "input_junctions",
        "skip_noncanonical_chrom",
        "skip_undefined_strand",
        "skip_low_unique_reads",
        "skip_annotated",
        "skip_no_protein_coding_overlap",
        "novel_candidates",
    ]
    details = " ".join(f"{field}={stats.get(field, 0)}" for field in fields)
    print(f"[spl2] {path}: {details}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    gtf = Path(args.gtf)
    star_root = Path(args.star_root)
    out_dir = Path(args.out_dir) if args.out_dir else None
    if not gtf.exists():
        raise SystemExit(f"ERROR: GTF does not exist: {gtf}")
    if not star_root.exists() or not star_root.is_dir():
        raise SystemExit(f"ERROR: STAR root does not exist or is not a directory: {star_root}")

    print(f"[spl2] loading GTF: {gtf}", file=sys.stderr)
    index = build_gtf_index(gtf, include_noncanonical=args.include_noncanonical)
    print(
        "[spl2] GTF index: "
        f"transcripts={index.transcript_count} "
        f"protein_coding_transcripts={index.protein_coding_transcript_count} "
        f"known_junctions={index.known_junction_count}",
        file=sys.stderr,
    )

    sj_files, duplicate_sj_files = discover_sj_files(star_root, args.sample_filter)
    if not sj_files:
        print(f"[spl2] no merged SJ.out.tab files found under {star_root}", file=sys.stderr)
        return 0
    for sample, ranked in sorted(duplicate_sj_files.items()):
        skipped = ", ".join(str(path) for path in ranked[1:])
        print(f"[spl2][dedupe] sample={sample} keeping {ranked[0]} skipped={skipped}", file=sys.stderr)

    written = 0
    skipped = 0
    for sj_path in sj_files:
        out_path = output_path_for(sj_path, star_root, out_dir, args.out_suffix, args.output_subdir)
        if out_path.exists() and not args.force:
            print(f"[spl2][skip] output exists: {out_path}", file=sys.stderr)
            skipped += 1
            continue
        print(f"[spl2] {out_path} <- {sj_path}", file=sys.stderr)
        if args.dry_run:
            written += 1
            continue
        stats = process_sj_file(
            sj_path=sj_path,
            out_path=out_path,
            index=index,
            min_unique_reads=args.min_unique_reads,
            include_noncanonical=args.include_noncanonical,
            keep_non_protein_coding=args.keep_non_protein_coding,
        )
        print_stats(out_path, stats)
        written += 1

    print(f"[spl2][done] written={written} skipped={skipped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
