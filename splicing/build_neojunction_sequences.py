#!/usr/bin/env python3
"""Build Arriba-like neosplicing sequence outputs from spl3 event annotations."""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Set, Tuple


INPUT_SUFFIX = ".spl3.event_annotated.tsv"
OUTPUT_SUFFIX = ".spl4.neojunctions.tsv"
NT_FASTA_SUFFIX = ".spl4.neojunctions.nt.fa"
AA_FASTA_SUFFIX = ".spl4.neojunctions.aa.fa"
CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class Feature:
    def __init__(self, start: int, end: int, lab: str, rank: int, phase: Optional[int] = None) -> None:
        self.start = start
        self.end = end
        self.lab = lab
        self.rank = rank
        self.phase = phase

    @property
    def idx(self) -> str:
        return f"{'E' if self.lab == 'exon' else 'I'}{self.rank}"


class Transcript:
    def __init__(
        self,
        transcript_id: str,
        gene_id: str,
        gene_name: str,
        gene_type: str,
        transcript_type: str,
        chrom: str,
        strand: str,
    ) -> None:
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.transcript_type = transcript_type
        self.chrom = chrom
        self.chrom_key = chrom_key(chrom)
        self.strand = strand
        self.exons_raw: List[Tuple[int, int, Optional[int]]] = []
        self.cds_raw: List[Tuple[int, int, Optional[int]]] = []
        self.features: List[Feature] = []
        self.cds_features: List[Feature] = []

    @property
    def is_protein_coding(self) -> bool:
        return self.gene_type == "protein_coding" or self.transcript_type == "protein_coding"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Build an Arriba-like neosplicing TSV and full nucleotide/amino-acid FASTA "
            "from spl3 SSNIP-style event annotations."
        )
    )
    ap.add_argument("--gtf", required=True, help="GENCODE/GTF annotation, optionally .gz")
    ap.add_argument("--fasta", required=True, help="Reference genome FASTA, optionally .gz")
    ap.add_argument("--root", required=True, help="Root to scan for spl3 *.event_annotated.tsv files")
    ap.add_argument(
        "--out-dir",
        default="",
        help="Optional root for spl4 outputs. Input paths relative to --root are preserved.",
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
    ap.add_argument("--sample-filter", action="append", default=[], help="Substring filter for input paths")
    ap.add_argument("--input-suffix", default=INPUT_SUFFIX, help="Suffix identifying spl3 input TSVs")
    ap.add_argument("--out-suffix", default=OUTPUT_SUFFIX, help="Suffix for the spl4 Arriba-like TSV")
    ap.add_argument("--nt-fasta-suffix", default=NT_FASTA_SUFFIX, help="Suffix for full nucleotide FASTA")
    ap.add_argument("--aa-fasta-suffix", default=AA_FASTA_SUFFIX, help="Suffix for full amino-acid FASTA")
    ap.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    ap.add_argument("--dry-run", action="store_true", help="Report planned outputs without writing")
    ap.add_argument("--include-noncanonical", action="store_true", help="Keep non-canonical contigs in the GTF")
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
    if value in {"", "."}:
        return None
    match = re.search(r"\d+", value)
    return int(match.group(0)) if match else None


def parse_phase(value: str) -> Optional[int]:
    if value in {"", "."}:
        return None
    try:
        phase = int(value)
    except ValueError:
        return None
    return phase if phase in {0, 1, 2} else None


def finalize_transcript(tx: Transcript) -> None:
    exons = sorted(tx.exons_raw, key=lambda item: (item[0], item[1]))
    raw_ranks = [rank for _, _, rank in exons]
    use_raw = all(rank is not None for rank in raw_ranks) and len(set(raw_ranks)) == len(raw_ranks)
    exon_ranks: Dict[Tuple[int, int], int] = {}
    if use_raw:
        for start, end, rank in exons:
            exon_ranks[(start, end)] = int(rank)  # type: ignore[arg-type]
    else:
        transcript_order = sorted(exons, key=lambda item: (item[0], item[1]), reverse=(tx.strand == "-"))
        for rank, (start, end, _) in enumerate(transcript_order, start=1):
            exon_ranks[(start, end)] = rank

    features: List[Feature] = []
    for start, end, _ in exons:
        features.append(Feature(start, end, "exon", exon_ranks[(start, end)]))
    for left, right in zip(exons, exons[1:]):
        intron_start = left[1] + 1
        intron_end = right[0] - 1
        if intron_start <= intron_end:
            features.append(Feature(intron_start, intron_end, "intron", min(exon_ranks[(left[0], left[1])], exon_ranks[(right[0], right[1])])))
    tx.features = sorted(features, key=lambda feature: (feature.start, feature.end, feature.lab))

    cds_features = []
    for start, end, phase in sorted(tx.cds_raw, key=lambda item: (item[0], item[1])):
        rank = exon_ranks.get(next(((e_start, e_end) for e_start, e_end, _ in exons if e_start <= start and end <= e_end), (start, end)), 0)
        cds_features.append(Feature(start, end, "cds", rank, phase))
    tx.cds_features = cds_features


def parse_transcript_ids_from_rows(rows: Iterable[Dict[str, str]]) -> Set[str]:
    ids: Set[str] = set()
    for row in rows:
        for key in ("matched_transcript_id", "transcript_id1", "transcript_ids"):
            value = row.get(key, "")
            if not value or value == "NA":
                continue
            for token in re.split(r"[;,|]", value):
                token = token.strip()
                if token and token != "NA":
                    ids.add(token)
                    ids.add(strip_version(token))
    return ids


def build_gtf_model(gtf: Path, wanted_transcripts: Set[str], include_noncanonical: bool) -> Dict[str, Transcript]:
    transcripts: Dict[str, Transcript] = {}

    with open_text(gtf) as fh:
        for lineno, raw_line in enumerate(fh, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue
            cols = raw_line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in {"exon", "CDS"}:
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
            if wanted_transcripts and transcript_id not in wanted_transcripts and strip_version(transcript_id) not in wanted_transcripts:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError as exc:
                raise SystemExit(f"ERROR: non-integer GTF coordinate at {gtf}:{lineno}: {exc}") from exc

            tx = transcripts.get(transcript_id)
            if tx is None:
                tx = Transcript(
                    transcript_id=transcript_id,
                    gene_id=attrs.get("gene_id", ""),
                    gene_name=attrs.get("gene_name", attrs.get("gene_id", "")),
                    gene_type=attrs.get("gene_type", attrs.get("gene_biotype", "")),
                    transcript_type=attrs.get("transcript_type", attrs.get("transcript_biotype", "")),
                    chrom=ssnip_chrom(chrom),
                    strand=strand,
                )
                transcripts[transcript_id] = tx

            if cols[2] == "exon":
                tx.exons_raw.append((start, end, parse_optional_int(attrs.get("exon_number", ""))))
            elif cols[2] == "CDS":
                tx.cds_raw.append((start, end, parse_phase(cols[7])))

    finalized: Dict[str, Transcript] = {}
    for tx in transcripts.values():
        if not tx.exons_raw or not tx.is_protein_coding:
            continue
        finalize_transcript(tx)
        finalized[tx.transcript_id] = tx
        finalized[strip_version(tx.transcript_id)] = tx
    return finalized


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

    # spl3 may be written directly under splicing_outdir as RNA_TUMOUR.spl3...
    # when earlier roots pointed at a sample folder, so match by TSV contents too.
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                haystack = " ".join(
                    row.get(column, "")
                    for column in (
                        "sample",
                        "source_sj",
                        "source_spl3",
                        "neojunction_id",
                        "junc_id",
                        "matched_gene_name",
                        "matched_transcript_id",
                    )
                ).lower()
                if any(filter_value in haystack for filter_value in lowered_filters):
                    return True
    except OSError:
        return False
    return False


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
    out_name = f"{input_path.name[:-len(input_suffix)]}{out_suffix}" if input_path.name.endswith(input_suffix) else f"{input_path.stem}{out_suffix}"
    if out_dir is None:
        return input_path.with_name(out_name)
    rel_parent = scoped_rel_parent(input_path.parent, root, output_subdir)
    return out_dir / rel_parent / out_name


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def collect_sequence_intervals(transcripts: Dict[str, Transcript]) -> DefaultDict[str, Set[Tuple[int, int]]]:
    intervals: DefaultDict[str, Set[Tuple[int, int]]] = defaultdict(set)
    seen: Set[str] = set()
    for tx in transcripts.values():
        if tx.transcript_id in seen:
            continue
        seen.add(tx.transcript_id)
        for feature in tx.features:
            intervals[tx.chrom].add((feature.start, feature.end))
        for feature in tx.cds_features:
            intervals[tx.chrom].add((feature.start, feature.end))
    return intervals


def fasta_aliases(chrom: str) -> Set[str]:
    aliases = {chrom, chrom_key(chrom)}
    if not chrom.startswith("chr"):
        aliases.add(f"chr{chrom}")
    if chrom_key(chrom) == "MT":
        aliases.update({"M", "MT", "chrM"})
    return aliases


def stream_fasta_intervals(fasta: Path, intervals_by_chrom: DefaultDict[str, Set[Tuple[int, int]]]) -> Dict[Tuple[str, int, int], str]:
    wanted_by_alias: Dict[str, Tuple[str, List[Tuple[int, int]]]] = {}
    for chrom, intervals in intervals_by_chrom.items():
        interval_list = sorted(intervals)
        for alias in fasta_aliases(chrom):
            wanted_by_alias[alias] = (chrom, interval_list)

    chunks: DefaultDict[Tuple[str, int, int], List[str]] = defaultdict(list)
    current_name = ""
    current_chrom = ""
    current_intervals: List[Tuple[int, int]] = []
    line_start = 1
    pointer = 0

    with open_text(fasta) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                current_chrom = ""
                current_intervals = []
                line_start = 1
                pointer = 0
                if current_name in wanted_by_alias:
                    current_chrom, current_intervals = wanted_by_alias[current_name]
                continue
            if not current_intervals:
                continue
            seq = line.strip()
            line_end = line_start + len(seq) - 1
            while pointer < len(current_intervals) and current_intervals[pointer][1] < line_start:
                pointer += 1
            idx = pointer
            while idx < len(current_intervals) and current_intervals[idx][0] <= line_end:
                start, end = current_intervals[idx]
                overlap_start = max(start, line_start)
                overlap_end = min(end, line_end)
                if overlap_start <= overlap_end:
                    offset_start = overlap_start - line_start
                    offset_end = overlap_end - line_start + 1
                    chunks[(current_chrom, start, end)].append(seq[offset_start:offset_end].upper())
                idx += 1
            line_start = line_end + 1

    seqs = {key: "".join(parts) for key, parts in chunks.items()}
    missing = []
    for chrom, intervals in intervals_by_chrom.items():
        for start, end in intervals:
            seq = seqs.get((chrom, start, end), "")
            if len(seq) != end - start + 1:
                missing.append(f"{chrom}:{start}-{end}")
                if len(missing) >= 5:
                    break
        if len(missing) >= 5:
            break
    if missing:
        raise SystemExit(f"ERROR: failed to fetch FASTA intervals, examples: {', '.join(missing)}")
    return seqs


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def interval_sequence(
    seq_cache: Dict[Tuple[str, int, int], str],
    tx: Transcript,
    feature: Feature,
    start: int,
    end: int,
) -> Tuple[str, List[int]]:
    full = seq_cache[(tx.chrom, feature.start, feature.end)]
    piece = full[start - feature.start : end - feature.start + 1]
    coords = list(range(start, end + 1))
    if tx.strand == "-":
        return reverse_complement(piece), list(reversed(coords))
    return piece.upper(), coords


def feature_contains(feature: Feature, pos: int) -> bool:
    return feature.start <= pos <= feature.end


def add_segment(
    seq_parts: List[str],
    coord_parts: List[int],
    seq_cache: Dict[Tuple[str, int, int], str],
    tx: Transcript,
    feature: Feature,
    start: int,
    end: int,
) -> None:
    if start > end:
        return
    seq, coords = interval_sequence(seq_cache, tx, feature, start, end)
    seq_parts.append(seq)
    coord_parts.extend(coords)


def transcript_ordered_features(tx: Transcript) -> List[Feature]:
    return sorted(tx.features, key=lambda feature: (feature.start, feature.end), reverse=(tx.strand == "-"))


def build_prefix(
    tx: Transcript,
    donor_pos: int,
    seq_cache: Dict[Tuple[str, int, int], str],
) -> Tuple[str, List[int]]:
    seq_parts: List[str] = []
    coords: List[int] = []
    for feature in transcript_ordered_features(tx):
        if tx.strand == "+":
            if feature.end < donor_pos:
                if feature.lab == "exon":
                    add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, feature.end)
                continue
            if feature_contains(feature, donor_pos):
                add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, donor_pos)
                break
            if feature.start > donor_pos:
                break
        else:
            if feature.start > donor_pos:
                if feature.lab == "exon":
                    add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, feature.end)
                continue
            if feature_contains(feature, donor_pos):
                add_segment(seq_parts, coords, seq_cache, tx, feature, donor_pos, feature.end)
                break
            if feature.end < donor_pos:
                break
    return "".join(seq_parts), coords


def build_suffix(
    tx: Transcript,
    acceptor_pos: int,
    seq_cache: Dict[Tuple[str, int, int], str],
) -> Tuple[str, List[int]]:
    seq_parts: List[str] = []
    coords: List[int] = []
    started = False
    for feature in transcript_ordered_features(tx):
        if tx.strand == "+":
            if feature.end < acceptor_pos:
                continue
            if feature_contains(feature, acceptor_pos):
                add_segment(seq_parts, coords, seq_cache, tx, feature, acceptor_pos, feature.end)
                started = True
                continue
            if started and feature.lab == "exon":
                add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, feature.end)
            elif feature.start > acceptor_pos and feature.lab == "exon":
                add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, feature.end)
                started = True
        else:
            if feature.start > acceptor_pos:
                continue
            if feature_contains(feature, acceptor_pos):
                add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, acceptor_pos)
                started = True
                continue
            if started and feature.lab == "exon":
                add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, feature.end)
            elif feature.end < acceptor_pos and feature.lab == "exon":
                add_segment(seq_parts, coords, seq_cache, tx, feature, feature.start, feature.end)
                started = True
    return "".join(seq_parts), coords


def parse_coords(row: Dict[str, str]) -> Tuple[str, str, int, int]:
    chrom = row.get("chrom", "")
    strand = row.get("strand", "")
    if row.get("star_intron_start") and row.get("star_intron_end"):
        return chrom, strand, int(row["star_intron_start"]), int(row["star_intron_end"])
    junc_id = row.get("junc_id", "")
    chrom_part, strand_part, coord_part = junc_id.split(":", 2)
    left, right = coord_part.split("-", 1)
    return chrom_part, strand_part, int(left) + 1, int(right)


def translate(seq: str) -> str:
    aa = []
    clean = seq.upper().replace("U", "T")
    for i in range(0, len(clean) - 2, 3):
        aa.append(CODON_TABLE.get(clean[i : i + 3], "X"))
    return "".join(aa)


def trim_at_stop(aa: str) -> Tuple[str, bool]:
    if "*" not in aa:
        return aa, False
    return aa.split("*", 1)[0], True


def longest_frame_translation(nt_seq: str) -> Tuple[str, int]:
    best_aa = ""
    best_frame = 0
    for frame in range(3):
        aa, _ = trim_at_stop(translate(nt_seq[frame:]))
        if len(aa) > len(best_aa):
            best_aa = aa
            best_frame = frame
    return best_aa, best_frame


def cds_anchor(tx: Transcript) -> Optional[int]:
    if not tx.cds_features:
        return None
    ordered = sorted(tx.cds_features, key=lambda feature: (feature.start, feature.end), reverse=(tx.strand == "-"))
    first = ordered[0]
    phase = first.phase or 0
    if tx.strand == "+":
        return first.start + phase
    return first.end - phase


def mark_with_pipe(seq: str, offset: int) -> str:
    offset = max(0, min(offset, len(seq)))
    return f"{seq[:offset]}|{seq[offset:]}"


def safe_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def make_event_id(row: Dict[str, str], gene_name: str, transcript_id: str) -> str:
    junc_id = row.get("junc_id") or f"{row.get('chrom', 'NA')}:{row.get('strand', 'NA')}:{row.get('star_intron_start', 'NA')}-{row.get('star_intron_end', 'NA')}"
    clean = re.sub(r"[^A-Za-z0-9_.:+-]", "_", junc_id)
    return f"{gene_name}|{gene_name}_{clean}|{transcript_id}"


def choose_transcript(row: Dict[str, str], transcripts: Dict[str, Transcript]) -> Optional[Transcript]:
    candidates = []
    for key in ("matched_transcript_id", "transcript_id1", "transcript_ids"):
        value = row.get(key, "")
        if not value or value == "NA":
            continue
        candidates.extend(token for token in re.split(r"[;,|]", value) if token and token != "NA")
    for transcript_id in candidates:
        tx = transcripts.get(transcript_id) or transcripts.get(strip_version(transcript_id))
        if tx is not None:
            return tx
    return None


def build_sequences(
    row: Dict[str, str],
    tx: Transcript,
    seq_cache: Dict[Tuple[str, int, int], str],
) -> Dict[str, str]:
    _, strand, intron_start, intron_end = parse_coords(row)
    if strand != tx.strand:
        return {"error": "strand_mismatch"}

    if tx.strand == "+":
        donor_pos = intron_start - 1
        acceptor_pos = intron_end + 1
    else:
        donor_pos = intron_end + 1
        acceptor_pos = intron_start - 1

    left_seq, left_coords = build_prefix(tx, donor_pos, seq_cache)
    right_seq, right_coords = build_suffix(tx, acceptor_pos, seq_cache)
    nt_seq = f"{left_seq}{right_seq}"
    coords = left_coords + right_coords
    junction_nt_offset = len(left_seq)
    nt_marked = mark_with_pipe(nt_seq, junction_nt_offset)

    flags = []
    anchor = cds_anchor(tx)
    if anchor is not None and anchor in coords:
        cds_offset = coords.index(anchor)
        aa_raw = translate(nt_seq[cds_offset:])
        aa_seq, saw_stop = trim_at_stop(aa_raw)
        translation_frame = str(cds_offset % 3)
        if saw_stop:
            flags.append("stop_codon_after_junction_or_cds")
    else:
        aa_seq, fallback_frame = longest_frame_translation(nt_seq)
        cds_offset = fallback_frame
        translation_frame = f"fallback_frame_{fallback_frame}"
        flags.append("no_cds_anchor")

    if not aa_seq:
        flags.append("no_aa_sequence")

    rel_junction = junction_nt_offset - cds_offset
    if rel_junction < 0:
        junction_aa_offset = 0
        flags.append("junction_before_translation_start")
    else:
        junction_aa_offset = rel_junction // 3
        if rel_junction % 3 != 0:
            flags.append("junction_inside_codon")
    aa_marked = mark_with_pipe(aa_seq, junction_aa_offset) if aa_seq else "NA"
    reading_frame = "unknown"
    if aa_seq:
        reading_frame = "out_of_frame" if "junction_inside_codon" in flags else "inframe"

    return {
        "nt_sequence": nt_seq,
        "nt_sequence_junction_marked": nt_marked,
        "aa_sequence": aa_seq or "NA",
        "aa_sequence_junction_marked": aa_marked,
        "junction_nt_offset": str(junction_nt_offset),
        "junction_aa_offset": str(junction_aa_offset),
        "translation_frame": translation_frame,
        "reading_frame": reading_frame,
        "donor_breakpoint": f"{tx.chrom}:{donor_pos}",
        "acceptor_breakpoint": f"{tx.chrom}:{acceptor_pos}",
        "qc_flags": ";".join(flags) if flags else "PASS",
    }


def confidence(row: Dict[str, str], sequence_fields: Dict[str, str]) -> Tuple[str, str]:
    flags = [] if sequence_fields.get("qc_flags") == "PASS" else sequence_fields.get("qc_flags", "").split(";")
    unique_reads = safe_int(row.get("unique_reads", "0"))
    total_reads = safe_int(row.get("total_reads", "0"))
    overhang = safe_int(row.get("max_splice_overhang", "0"))
    event_class = row.get("ssnip_event_class", row.get("event_class", ""))

    if unique_reads < 10:
        flags.append("low_unique_reads")
    if total_reads < 10:
        flags.append("low_total_reads")
    if overhang < 20:
        flags.append("low_splice_overhang")
    if event_class in {"other", "no_transcript_match", ""}:
        flags.append("weak_event_class")
    if sequence_fields.get("aa_sequence_junction_marked", "NA") == "NA":
        flags.append("no_pipe_marked_aa")

    severe = {"no_aa_sequence", "strand_mismatch", "no_transcript_model", "no_pipe_marked_aa"}
    if any(flag in severe for flag in flags):
        level = "low"
    elif unique_reads >= 10 and total_reads >= 20 and overhang >= 20 and event_class not in {"other", "no_transcript_match", ""}:
        level = "high"
    elif unique_reads >= 10 and total_reads >= 10:
        level = "medium"
    else:
        level = "low"
    return level, ";".join(sorted(set(flag for flag in flags if flag))) if flags else "PASS"


def fasta_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:+|-]", "_", value)


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def process_file(
    input_path: Path,
    out_tsv: Path,
    out_nt_fasta: Path,
    out_aa_fasta: Path,
    transcripts: Dict[str, Transcript],
    seq_cache: Dict[Tuple[str, int, int], str],
) -> Dict[str, int]:
    rows = read_tsv(input_path)
    out_rows: List[Dict[str, str]] = []
    nt_records: List[Tuple[str, str]] = []
    aa_records: List[Tuple[str, str]] = []
    stats: DefaultDict[str, int] = defaultdict(int)

    for idx, row in enumerate(rows, start=1):
        stats["input_rows"] += 1
        tx = choose_transcript(row, transcripts)
        if tx is None:
            sequence_fields = {
                "nt_sequence": "NA",
                "nt_sequence_junction_marked": "NA",
                "aa_sequence": "NA",
                "aa_sequence_junction_marked": "NA",
                "junction_nt_offset": "NA",
                "junction_aa_offset": "NA",
                "translation_frame": "NA",
                "reading_frame": "unknown",
                "donor_breakpoint": "NA",
                "acceptor_breakpoint": "NA",
                "qc_flags": "no_transcript_model",
            }
            gene_id = row.get("matched_gene_id", row.get("gene_ids", "NA"))
            gene_name = row.get("matched_gene_name", row.get("gene_names", "NA"))
            transcript_id = row.get("matched_transcript_id", "NA")
        else:
            sequence_fields = build_sequences(row, tx, seq_cache)
            if "error" in sequence_fields:
                sequence_fields.update({
                    "nt_sequence": "NA",
                    "nt_sequence_junction_marked": "NA",
                    "aa_sequence": "NA",
                    "aa_sequence_junction_marked": "NA",
                    "junction_nt_offset": "NA",
                    "junction_aa_offset": "NA",
                    "translation_frame": "NA",
                    "reading_frame": "unknown",
                    "donor_breakpoint": "NA",
                    "acceptor_breakpoint": "NA",
                    "qc_flags": sequence_fields["error"],
                })
            gene_id = tx.gene_id or row.get("matched_gene_id", "NA")
            gene_name = tx.gene_name or row.get("matched_gene_name", "NA")
            transcript_id = tx.transcript_id

        conf, qc_flags = confidence(row, sequence_fields)
        sequence_fields["qc_flags"] = qc_flags
        chrom, strand, intron_start, intron_end = parse_coords(row)
        left_breakpoint = f"{chrom}:{intron_start - 1}"
        right_breakpoint = f"{chrom}:{intron_end + 1}"
        event_id = make_event_id(row, gene_name, transcript_id)
        nt_id = fasta_id(f"{event_id}|nt|row{idx}")
        aa_id = fasta_id(f"{event_id}|aa|row{idx}")

        if sequence_fields["nt_sequence"] != "NA":
            nt_records.append((nt_id, sequence_fields["nt_sequence"]))
        if sequence_fields["aa_sequence"] != "NA":
            aa_records.append((aa_id, sequence_fields["aa_sequence"]))

        unique_reads = row.get("unique_reads", "0") or "0"
        multimap_reads = row.get("multimap_reads", "0") or "0"
        total_reads = row.get("total_reads", unique_reads) or unique_reads
        event_class = row.get("ssnip_event_class", row.get("event_class", "NA"))
        peptide_sequence = sequence_fields["aa_sequence_junction_marked"]

        out = {
            "#gene1": gene_name,
            "gene2": gene_name,
            "gene_id1": gene_id,
            "gene_id2": gene_id,
            "transcript_id1": transcript_id,
            "transcript_id2": transcript_id,
            "discordant_mates": unique_reads,
            "split_reads1": unique_reads,
            "split_reads2": "0",
            "reading_frame": sequence_fields["reading_frame"],
            "site1": row.get("left_feature", "NA"),
            "site2": row.get("right_feature", "NA"),
            "breakpoint1": left_breakpoint,
            "breakpoint2": right_breakpoint,
            "peptide_sequence": peptide_sequence,
            "confidence": conf,
            "sample": row.get("sample", "NA"),
            "event_type": "NEOSPLICING",
            "source_set": "NEOSPLICING",
            "event_id": event_id,
            "neojunction_id": row.get("junc_id", "NA"),
            "fusion_id": event_id,
            "chrom": chrom,
            "strand": strand,
            "star_intron_start": str(intron_start),
            "star_intron_end": str(intron_end),
            "donor_boundary": row.get("donor_boundary", sequence_fields.get("donor_breakpoint", "NA").split(":")[-1]),
            "acceptor_boundary": row.get("acceptor_boundary", sequence_fields.get("acceptor_breakpoint", "NA").split(":")[-1]),
            "ssnip_event_class": event_class,
            "ssnip_event_raw": row.get("ssnip_event_raw", "NA"),
            "left_feature": row.get("left_feature", "NA"),
            "right_feature": row.get("right_feature", "NA"),
            "skipped_exon_ids": row.get("skipped_exon_ids", "NA"),
            "unique_reads": unique_reads,
            "multimap_reads": multimap_reads,
            "total_reads": total_reads,
            "max_splice_overhang": row.get("max_splice_overhang", "NA"),
            "qc_confidence": conf,
            "qc_flags": qc_flags,
            "nt_sequence_junction_marked": sequence_fields["nt_sequence_junction_marked"],
            "aa_sequence_junction_marked": sequence_fields["aa_sequence_junction_marked"],
            "junction_nt_offset": sequence_fields["junction_nt_offset"],
            "junction_aa_offset": sequence_fields["junction_aa_offset"],
            "translation_frame": sequence_fields["translation_frame"],
            "nt_fasta_id": nt_id if sequence_fields["nt_sequence"] != "NA" else "NA",
            "aa_fasta_id": aa_id if sequence_fields["aa_sequence"] != "NA" else "NA",
            "source_sj": row.get("source_sj", "NA"),
            "source_spl3": str(input_path),
        }
        out_rows.append(out)
        stats[f"confidence_{conf}"] += 1

    fieldnames = [
        "#gene1", "gene2", "gene_id1", "gene_id2", "transcript_id1", "transcript_id2",
        "discordant_mates", "split_reads1", "split_reads2", "reading_frame", "site1", "site2",
        "breakpoint1", "breakpoint2", "peptide_sequence", "confidence",
        "sample", "event_type", "source_set", "event_id", "neojunction_id", "fusion_id",
        "chrom", "strand", "star_intron_start", "star_intron_end", "donor_boundary", "acceptor_boundary",
        "ssnip_event_class", "ssnip_event_raw", "left_feature", "right_feature", "skipped_exon_ids",
        "unique_reads", "multimap_reads", "total_reads", "max_splice_overhang",
        "qc_confidence", "qc_flags", "nt_sequence_junction_marked", "aa_sequence_junction_marked",
        "junction_nt_offset", "junction_aa_offset", "translation_frame", "nt_fasta_id", "aa_fasta_id",
        "source_sj", "source_spl3",
    ]

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)
    with out_nt_fasta.open("w", encoding="utf-8") as out:
        for rec_id, seq in nt_records:
            out.write(f">{rec_id}\n{wrap_fasta(seq)}\n")
    with out_aa_fasta.open("w", encoding="utf-8") as out:
        for rec_id, seq in aa_records:
            out.write(f">{rec_id}\n{wrap_fasta(seq)}\n")
    stats["output_rows"] = len(out_rows)
    stats["nt_fasta_records"] = len(nt_records)
    stats["aa_fasta_records"] = len(aa_records)
    return dict(stats)


def output_paths(input_path: Path, root: Path, out_dir: Optional[Path], args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    out_tsv = output_path_for(input_path, root, out_dir, args.input_suffix, args.out_suffix, args.output_subdir)
    out_nt = output_path_for(input_path, root, out_dir, args.input_suffix, args.nt_fasta_suffix, args.output_subdir)
    out_aa = output_path_for(input_path, root, out_dir, args.input_suffix, args.aa_fasta_suffix, args.output_subdir)
    return out_tsv, out_nt, out_aa


def print_stats(path: Path, stats: Dict[str, int]) -> None:
    fields = ["input_rows", "output_rows", "nt_fasta_records", "aa_fasta_records", "confidence_high", "confidence_medium", "confidence_low"]
    details = " ".join(f"{field}={stats.get(field, 0)}" for field in fields)
    print(f"[spl4] {path}: {details}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    gtf = Path(args.gtf)
    fasta = Path(args.fasta)
    root = Path(args.root)
    out_dir = Path(args.out_dir) if args.out_dir else None
    if not gtf.exists():
        raise SystemExit(f"ERROR: GTF does not exist: {gtf}")
    if not fasta.exists():
        raise SystemExit(f"ERROR: FASTA does not exist: {fasta}")
    if not root.exists() or not root.is_dir():
        raise SystemExit(f"ERROR: spl3 root does not exist or is not a directory: {root}")

    inputs = discover_inputs(root, args.sample_filter, args.input_suffix)
    if not inputs:
        print(f"[spl4] no spl3 TSV files found under {root}", file=sys.stderr)
        return 0

    cached_rows = {path: read_tsv(path) for path in inputs}
    wanted_transcripts: Set[str] = set()
    for rows in cached_rows.values():
        wanted_transcripts.update(parse_transcript_ids_from_rows(rows))

    print(f"[spl4] loading GTF: {gtf}", file=sys.stderr)
    transcripts = build_gtf_model(gtf, wanted_transcripts, include_noncanonical=args.include_noncanonical)
    print(f"[spl4] GTF protein-coding transcript models={len({tx.transcript_id for tx in transcripts.values()})}", file=sys.stderr)

    intervals = collect_sequence_intervals(transcripts)
    print(f"[spl4] loading FASTA intervals from: {fasta}", file=sys.stderr)
    seq_cache = stream_fasta_intervals(fasta, intervals)
    print(f"[spl4] FASTA intervals loaded={len(seq_cache)}", file=sys.stderr)

    written = 0
    skipped = 0
    for input_path in inputs:
        out_tsv, out_nt, out_aa = output_paths(input_path, root, out_dir, args)
        outputs_exist = out_tsv.exists() and out_nt.exists() and out_aa.exists()
        if outputs_exist and not args.force:
            print(f"[spl4][skip] outputs exist: {out_tsv}", file=sys.stderr)
            skipped += 1
            continue
        print(f"[spl4] {out_tsv} <- {input_path}", file=sys.stderr)
        if args.dry_run:
            written += 1
            continue
        stats = process_file(input_path, out_tsv, out_nt, out_aa, transcripts, seq_cache)
        print_stats(out_tsv, stats)
        written += 1

    print(f"[spl4][done] written={written} skipped={skipped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
