#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, TextIO, Tuple


PRIMARY_GTF_FEATURES = {
    "CDS",
    "UTR",
    "exon",
    "five_prime_utr",
    "start_codon",
    "stop_codon",
    "three_prime_utr",
}
FALLBACK_GTF_FEATURES = {"transcript"}


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom(chrom: str) -> str:
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def open_alignment(path: str):
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            "ERROR: pysam is required for RNA strandedness filtering. "
            "Load a Python environment that includes pysam in modules_rna/research_python_modules."
        ) from exc
    return pysam.AlignmentFile(path, "rb")


def build_contig_lookup(references: Sequence[str]) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for ref in references:
        norm = normalize_chrom(ref)
        lookup.setdefault(norm, ref)
        lookup.setdefault(ref, ref)
    return lookup


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


def format_info(original_info: str, extra_fields: Sequence[Tuple[str, str]]) -> str:
    parts: List[str] = []
    if original_info not in ("", "."):
        parts.extend(original_info.split(";"))
    parts.extend(f"{key}={value}" for key, value in extra_fields if value not in ("", ".", "NA"))
    return ";".join(parts) if parts else "."


def protocol_alias(protocol: str) -> str:
    text = str(protocol).strip().lower()
    aliases = {
        "fr-firststrand": "fr-firststrand",
        "firststrand": "fr-firststrand",
        "reverse": "fr-firststrand",
        "reverse-stranded": "fr-firststrand",
        "1+-,1-+,2++,2--": "fr-firststrand",
        "fr-secondstrand": "fr-secondstrand",
        "secondstrand": "fr-secondstrand",
        "forward": "fr-secondstrand",
        "forward-stranded": "fr-secondstrand",
        "1++,1--,2+-,2-+": "fr-secondstrand",
    }
    if text not in aliases:
        raise SystemExit(f"ERROR: unsupported protocol '{protocol}'")
    return aliases[text]


def transcript_strand_from_read(aln, protocol: str) -> str:
    protocol = protocol_alias(protocol)
    if aln.is_paired:
        if aln.is_read1:
            if protocol == "fr-firststrand":
                return "+" if aln.is_reverse else "-"
            return "-" if aln.is_reverse else "+"
        if aln.is_read2:
            if protocol == "fr-firststrand":
                return "-" if aln.is_reverse else "+"
            return "+" if aln.is_reverse else "-"
    if protocol == "fr-firststrand":
        return "+" if aln.is_reverse else "-"
    return "-" if aln.is_reverse else "+"


def maybe_num(value: Optional[float], digits: int = 6) -> str:
    if value is None:
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def join_values(values: Iterable[str]) -> str:
    clean = sorted({str(v) for v in values if str(v) not in ("", ".", "NA")})
    return "|".join(clean) if clean else "NA"


def add_unique(store: Dict[str, Set[str]], key: str, value: str) -> None:
    text = str(value).strip()
    if text and text not in (".", "-"):
        store[key].add(text)


def classify_edit_sig(info_map: Dict[str, str]) -> str:
    labels: List[str] = []
    edit_sig = info_map.get("EDIT_SIG", "")
    for token in str(edit_sig).split(","):
        token = token.strip()
        if token:
            labels.append(token)
    if "ADAR_SIG" in info_map and "ADAR" not in labels:
        labels.append("ADAR")
    if "APOBEC3_SIG" in info_map and "APOBEC3" not in labels:
        labels.append("APOBEC3")
    if not labels:
        return "NA"
    return "|".join(sorted(set(labels)))


def load_gtf_bins(gtf_path: str, bin_size: int = 100_000):
    primary_index: Dict[str, Dict[int, List[Tuple[int, int, str]]]] = defaultdict(lambda: defaultdict(list))
    fallback_index: Dict[str, Dict[int, List[Tuple[int, int, str]]]] = defaultdict(lambda: defaultdict(list))

    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = normalize_chrom(fields[0])
            feature = fields[2]
            if feature in PRIMARY_GTF_FEATURES:
                target = primary_index
            elif feature in FALLBACK_GTF_FEATURES:
                target = fallback_index
            else:
                continue
            strand = fields[6]
            if strand not in ("+", "-"):
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            first_bin = start // bin_size
            last_bin = (end - 1) // bin_size
            for bin_id in range(first_bin, last_bin + 1):
                target[chrom][bin_id].append((start, end, strand))

    return {"primary": primary_index, "fallback": fallback_index, "bin_size": bin_size}


def query_transcript_strands(gtf_bins, chrom: str, pos: int) -> Tuple[Set[str], str]:
    norm_chrom = normalize_chrom(chrom)
    pos0 = pos - 1
    bin_id = pos0 // gtf_bins["bin_size"]

    def collect(index: Dict[str, Dict[int, List[Tuple[int, int, str]]]]) -> Set[str]:
        strands: Set[str] = set()
        for start, end, strand in index.get(norm_chrom, {}).get(bin_id, []):
            if start <= pos0 < end:
                strands.add(strand)
        return strands

    strands = collect(gtf_bins["primary"])
    if strands:
        return strands, "exonic"
    strands = collect(gtf_bins["fallback"])
    if strands:
        return strands, "transcript"
    return set(), "missing"


def count_alt_reads(
    bam,
    actual_chrom: str,
    pos: int,
    alt: str,
    expected_strand: str,
    protocol: str,
    min_mapq: int,
    min_baseq: int,
) -> Tuple[int, int, int, Optional[float]]:
    expected = 0
    opposite = 0
    total = 0

    for pileup_col in bam.pileup(
        actual_chrom,
        pos - 1,
        pos,
        truncate=True,
        stepper="samtools",
        min_base_quality=min_baseq,
        min_mapping_quality=min_mapq,
    ):
        if pileup_col.reference_pos != pos - 1:
            continue
        for pileup_read in pileup_col.pileups:
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            aln = pileup_read.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.is_duplicate or aln.is_qcfail:
                continue
            qpos = pileup_read.query_position
            if qpos is None:
                continue
            base = aln.query_sequence[qpos]
            if base is None or base.upper() != alt.upper():
                continue
            if aln.query_qualities is not None:
                if qpos >= len(aln.query_qualities) or aln.query_qualities[qpos] < min_baseq:
                    continue
            tx_strand = transcript_strand_from_read(aln, protocol)
            total += 1
            if tx_strand == expected_strand:
                expected += 1
            else:
                opposite += 1
        break

    frac = (expected / total) if total > 0 else None
    return expected, opposite, total, frac


def derive_reason(
    transcript_strands: Set[str],
    total_alt: int,
    expected_frac: Optional[float],
    min_expected_frac: float,
    missing_bam_contig: bool = False,
) -> str:
    if not transcript_strands:
        return "MISSING_TRANSCRIPT_STRAND"
    if len(transcript_strands) > 1:
        return "AMBIGUOUS_TRANSCRIPT_STRAND"
    if missing_bam_contig:
        return "MISSING_BAM_CONTIG"
    if total_alt == 0:
        return "NO_ALT_READ_SUPPORT"
    if expected_frac is None:
        return "NO_EXPECTED_FRACTION"
    if expected_frac < min_expected_frac:
        return "OPPOSITE_STRAND_ENRICHED"
    return "PASS"


def write_tsv(path: str, rows: Sequence[Dict[str, str]], fieldnames: Sequence[str]) -> None:
    if not path:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    ap = argparse.ArgumentParser(description="Filter RNA-edit VCF by strandedness using RNA BAM and GTF transcript strands.")
    ap.add_argument("-i", "--input", required=True, help="Input RNA-edit VCF (preferably split, e.g. rna5 output)")
    ap.add_argument("-o", "--output", required=True, help="Output filtered VCF")
    ap.add_argument("--bam", required=True, help="RNA BAM")
    ap.add_argument("--gtf", required=True, help="Reference GTF/GTF.gz used for transcript strand lookup")
    ap.add_argument("--patient", default="NA", help="Patient/sample identifier for TSV sidecars")
    ap.add_argument("--protocol", default="fr-firststrand")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-expected-frac", type=float, default=0.8)
    ap.add_argument("--support-tsv", default="")
    ap.add_argument("--blacklist-tsv", default="")
    ap.add_argument("--stats", default="")
    ap.add_argument(
        "--filter-source-set",
        default="",
        help="Optional SOURCE_SET token to restrict filtering to; non-matching records are preserved unchanged",
    )
    args = ap.parse_args()

    protocol = protocol_alias(args.protocol)
    gtf_bins = load_gtf_bins(args.gtf)
    support_rows: List[Dict[str, str]] = []
    blacklist_rows: List[Dict[str, str]] = []
    reason_counts: Counter[str] = Counter()
    input_records = 0
    kept_records = 0
    preserved_nonmatching_records = 0

    with open_alignment(args.bam) as bam:
        contig_lookup = build_contig_lookup(bam.references)
        with open_text(args.input) as fin, open(args.output, "w") as fout:
            for line in fin:
                if line.startswith("##"):
                    fout.write(line)
                    continue
                if line.startswith("#CHROM"):
                    fout.write('##INFO=<ID=RNA_STRAND,Number=1,Type=String,Description="Expected transcript strand used for RNA strandedness filter (+ or -)">\n')
                    fout.write('##INFO=<ID=RNA_STRAND_SET,Number=1,Type=String,Description="Unique overlapping transcript strands from GTF at this position">\n')
                    fout.write('##INFO=<ID=RNA_STRAND_SOURCE,Number=1,Type=String,Description="Annotation tier providing RNA_STRAND (exonic, transcript, missing)">\n')
                    fout.write('##INFO=<ID=RNA_EXPECTED_ALT_READS,Number=1,Type=Integer,Description="ALT-supporting RNA reads consistent with expected transcript strand">\n')
                    fout.write('##INFO=<ID=RNA_OPPOSITE_ALT_READS,Number=1,Type=Integer,Description="ALT-supporting RNA reads on the opposite transcript strand">\n')
                    fout.write('##INFO=<ID=RNA_TOTAL_ALT_READS,Number=1,Type=Integer,Description="Total ALT-supporting RNA reads counted for strandedness QC">\n')
                    fout.write('##INFO=<ID=RNA_EXPECTED_ALT_FRAC,Number=1,Type=Float,Description="Fraction of ALT-supporting RNA reads on the expected transcript strand">\n')
                    fout.write('##INFO=<ID=RNA_STRAND_REASON,Number=1,Type=String,Description="Result of RNA strandedness QC (PASS or blacklist reason)">\n')
                    fout.write(line)
                    continue
                if not line or line.startswith("#"):
                    fout.write(line)
                    continue

                input_records += 1
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                chrom = normalize_chrom(cols[0])
                pos = int(cols[1])
                ref = cols[3]
                alt = cols[4]
                info_map = parse_info(cols[7])
                source_tokens = {token.strip() for token in str(info_map.get("SOURCE_SET", "")).split(",") if token.strip()}
                if args.filter_source_set and args.filter_source_set not in source_tokens:
                    fout.write(line)
                    kept_records += 1
                    preserved_nonmatching_records += 1
                    continue
                filter_value = cols[6] if cols[6] else "NA"
                edit_sig = classify_edit_sig(info_map)
                known_db = info_map.get("KNOWN_RNAEDIT_DB", "NA")

                transcript_strands, annotation_source = query_transcript_strands(gtf_bins, chrom, pos)
                expected_strand = sorted(transcript_strands)[0] if len(transcript_strands) == 1 else "NA"

                expected_alt = 0
                opposite_alt = 0
                total_alt = 0
                expected_frac: Optional[float] = None
                missing_bam_contig = False
                if len(transcript_strands) == 1:
                    actual_chrom = contig_lookup.get(chrom)
                    if not actual_chrom:
                        missing_bam_contig = True
                    else:
                        expected_alt, opposite_alt, total_alt, expected_frac = count_alt_reads(
                            bam=bam,
                            actual_chrom=actual_chrom,
                            pos=pos,
                            alt=alt,
                            expected_strand=expected_strand,
                            protocol=protocol,
                            min_mapq=args.min_mapq,
                            min_baseq=args.min_baseq,
                        )

                reason = derive_reason(
                    transcript_strands=transcript_strands,
                    total_alt=total_alt,
                    expected_frac=expected_frac,
                    min_expected_frac=args.min_expected_frac,
                    missing_bam_contig=missing_bam_contig,
                )
                blacklisted = "1" if reason != "PASS" else "0"
                reason_counts[reason] += 1

                support_row = {
                    "SAMPLE": args.patient,
                    "chrom": chrom,
                    "pos": str(pos),
                    "ref": ref,
                    "alt": alt,
                    "filter": filter_value,
                    "edit_sig": edit_sig,
                    "known_rnaedit_db": known_db if known_db else "NA",
                    "transcript_strand": expected_strand,
                    "transcript_strand_set": join_values(transcript_strands),
                    "strand_annotation_source": annotation_source,
                    "expected_alt_reads": maybe_num(expected_alt, digits=3),
                    "opposite_alt_reads": maybe_num(opposite_alt, digits=3),
                    "total_alt_reads": maybe_num(total_alt, digits=3),
                    "expected_alt_fraction": maybe_num(expected_frac, digits=6),
                    "strand_filter_reason": reason,
                    "blacklisted": blacklisted,
                }
                support_rows.append(support_row)
                if blacklisted == "1":
                    blacklist_rows.append(support_row)
                    continue

                extra_info = [
                    ("RNA_STRAND", expected_strand),
                    ("RNA_STRAND_SET", join_values(transcript_strands)),
                    ("RNA_STRAND_SOURCE", annotation_source),
                    ("RNA_EXPECTED_ALT_READS", maybe_num(expected_alt, digits=3)),
                    ("RNA_OPPOSITE_ALT_READS", maybe_num(opposite_alt, digits=3)),
                    ("RNA_TOTAL_ALT_READS", maybe_num(total_alt, digits=3)),
                    ("RNA_EXPECTED_ALT_FRAC", maybe_num(expected_frac, digits=6)),
                    ("RNA_STRAND_REASON", reason),
                ]
                cols[7] = format_info(cols[7], extra_info)
                fout.write("\t".join(cols) + "\n")
                kept_records += 1

    fieldnames = [
        "SAMPLE",
        "chrom",
        "pos",
        "ref",
        "alt",
        "filter",
        "edit_sig",
        "known_rnaedit_db",
        "transcript_strand",
        "transcript_strand_set",
        "strand_annotation_source",
        "expected_alt_reads",
        "opposite_alt_reads",
        "total_alt_reads",
        "expected_alt_fraction",
        "strand_filter_reason",
        "blacklisted",
    ]
    write_tsv(args.support_tsv, support_rows, fieldnames)
    write_tsv(args.blacklist_tsv, blacklist_rows, fieldnames)

    if args.stats:
        os.makedirs(os.path.dirname(args.stats) or ".", exist_ok=True)
        with open(args.stats, "w") as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"input_records\t{input_records}\n")
            fh.write(f"kept_records\t{kept_records}\n")
            fh.write(f"preserved_nonmatching_records\t{preserved_nonmatching_records}\n")
            fh.write(f"blacklisted_records\t{len(blacklist_rows)}\n")
            for reason, count in sorted(reason_counts.items()):
                fh.write(f"reason_{reason}\t{count}\n")


if __name__ == "__main__":
    main()
