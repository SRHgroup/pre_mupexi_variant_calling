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


def normalize_chrom(chrom):
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def open_alignment(path):
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            "ERROR: pysam is required for strand blacklist generation. "
            "Load a Python environment that includes pysam in research_python_modules."
        ) from exc
    return pysam.AlignmentFile(path, "rb")


def parse_info(info):
    out = {}
    if info in ("", "."):
        return out
    for token in info.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
        else:
            out[token] = True
    return out


def parse_vep_location(loc):
    text = str(loc).strip()
    if ":" not in text:
        return None, None
    chrom, rest = text.split(":", 1)
    pos_s = rest.split("-", 1)[0]
    try:
        pos = int(pos_s)
    except Exception:
        return None, None
    return normalize_chrom(chrom), pos


def parse_symbol(extra):
    m = re.search(r"(?:^|;)SYMBOL=([^;]+)", str(extra))
    return m.group(1) if m else ""


def add_unique(store, key, value):
    text = str(value).strip()
    if text and text not in (".", "-"):
        store[key].add(text)


def join_values(values):
    clean = [str(v) for v in values if str(v) not in ("", ".", "NA")]
    if not clean:
        return "NA"
    return "|".join(sorted(set(clean)))


def maybe_num(value, digits=6):
    if value is None:
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def build_contig_lookup(references):
    lookup = {}
    for ref in references:
        norm = normalize_chrom(ref)
        lookup.setdefault(norm, ref)
        lookup.setdefault(ref, ref)
    return lookup


def protocol_alias(protocol):
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


def transcript_strand_from_read(aln, protocol):
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
    # single-end fallback
    if protocol == "fr-firststrand":
        return "+" if aln.is_reverse else "-"
    return "-" if aln.is_reverse else "+"


def load_vep_map(path):
    out = {}
    if not path or not os.path.exists(path):
        return out

    header = None
    with open_text(path) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#Uploaded_variation"):
                header = [c.lstrip("#") for c in cols]
                continue
            if line.startswith("#"):
                continue

            if header:
                row = {header[i]: cols[i] if i < len(cols) else "" for i in range(len(header))}
            else:
                default_header = [
                    "Uploaded_variation",
                    "Location",
                    "Allele",
                    "Gene",
                    "Feature",
                    "Feature_type",
                    "Consequence",
                    "cDNA_position",
                    "CDS_position",
                    "Protein_position",
                    "Amino_acids",
                    "Codons",
                    "Existing_variation",
                    "Extra",
                ]
                row = {default_header[i]: cols[i] if i < len(cols) else "" for i in range(len(default_header))}

            chrom, pos = parse_vep_location(row.get("Location", ""))
            alt = str(row.get("Allele", "")).strip()
            if chrom is None or pos is None or not alt:
                continue

            key = (chrom, pos, alt)
            if key not in out:
                out[key] = defaultdict(set)

            rec = out[key]
            add_unique(rec, "uploaded_variation", row.get("Uploaded_variation", ""))
            add_unique(rec, "location", row.get("Location", ""))
            add_unique(rec, "gene_id", row.get("Gene", ""))
            add_unique(rec, "gene_symbol", parse_symbol(row.get("Extra", "")))
            add_unique(rec, "feature", row.get("Feature", ""))
            add_unique(rec, "feature_type", row.get("Feature_type", ""))
            add_unique(rec, "consequence", row.get("Consequence", ""))
            add_unique(rec, "existing_variation", row.get("Existing_variation", ""))
            strand_match = re.search(r"(?:^|;)STRAND=([+-]?[0-9]+)", str(row.get("Extra", "")))
            if strand_match:
                strand_val = strand_match.group(1)
                if strand_val == "1":
                    rec["transcript_strand"].add("+")
                elif strand_val == "-1":
                    rec["transcript_strand"].add("-")

    return out


def count_alt_reads(bam, actual_chrom, pos, alt, expected_strand, protocol, min_mapq, min_baseq):
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

    frac = None
    if total > 0:
        frac = expected / total
    return expected, opposite, total, frac


def iter_rna_edit_variants(vcf_path):
    with open_text(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")
            flt = fields[6] if fields[6] else "NA"
            info_map = parse_info(fields[7])
            source_set = str(info_map.get("SOURCE_SET", ""))
            if "RNA_EDIT" not in source_set:
                continue
            for alt in alts:
                yield {
                    "chrom": normalize_chrom(chrom),
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "filter": flt,
                    "source_set": source_set,
                }


def derive_reason(transcript_strands, total_alt, expected_alt, expected_frac, min_expected_frac, missing_bam_contig=False):
    if not transcript_strands:
        return "MISSING_TRANSCRIPT_STRAND", 1
    if len(transcript_strands) > 1:
        return "AMBIGUOUS_TRANSCRIPT_STRAND", 1
    if missing_bam_contig:
        return "MISSING_BAM_CONTIG", 1
    if total_alt == 0:
        return "NO_ALT_READ_SUPPORT", 1
    if expected_frac is None:
        return "NO_EXPECTED_FRACTION", 1
    if expected_frac < min_expected_frac:
        return "OPPOSITE_STRAND_ENRICHED", 1
    return "PASS", 0


def analyze_patient(patient, vcf_path, bam_path, vep_path, protocol, min_mapq, min_baseq, min_expected_frac):
    vep_map = load_vep_map(vep_path)
    rows = []
    with open_alignment(bam_path) as bam:
        contig_lookup = build_contig_lookup(bam.references)

        for var in iter_rna_edit_variants(vcf_path):
            key = (var["chrom"], var["pos"], var["alt"])
            ann = vep_map.get(key, {})
            strands = sorted(ann.get("transcript_strand", set()))
            expected_strand = strands[0] if len(strands) == 1 else "NA"

            expected_alt = 0
            opposite_alt = 0
            total_alt = 0
            expected_frac = None
            missing_bam_contig = False
            if len(strands) == 1:
                actual_chrom = contig_lookup.get(var["chrom"])
                if not actual_chrom:
                    missing_bam_contig = True
                else:
                    expected_alt, opposite_alt, total_alt, expected_frac = count_alt_reads(
                        bam=bam,
                        actual_chrom=actual_chrom,
                        pos=var["pos"],
                        alt=var["alt"],
                        expected_strand=expected_strand,
                        protocol=protocol,
                        min_mapq=min_mapq,
                        min_baseq=min_baseq,
                    )

            reason, blacklisted = derive_reason(
                strands,
                total_alt,
                expected_alt,
                expected_frac,
                min_expected_frac,
                missing_bam_contig=missing_bam_contig,
            )
            rows.append({
                "SAMPLE": patient,
                "uploaded_variation": join_values(ann.get("uploaded_variation", {f"{var['chrom']}_{var['pos']}_{var['ref']}/{var['alt']}"})),
                "location": join_values(ann.get("location", {f"{var['chrom']}:{var['pos']}"})),
                "chrom": var["chrom"],
                "pos": str(var["pos"]),
                "ref": var["ref"],
                "alt": var["alt"],
                "gene_id": join_values(ann.get("gene_id", set())),
                "gene_symbol": join_values(ann.get("gene_symbol", set())),
                "feature": join_values(ann.get("feature", set())),
                "feature_type": join_values(ann.get("feature_type", set())),
                "consequence": join_values(ann.get("consequence", set())),
                "existing_variation": join_values(ann.get("existing_variation", set())),
                "source_set": var["source_set"],
                "filter": var["filter"],
                "protocol": protocol_alias(protocol),
                "transcript_strand": expected_strand,
                "transcript_strand_set": join_values(strands),
                "expected_alt_reads": maybe_num(expected_alt, digits=3),
                "opposite_alt_reads": maybe_num(opposite_alt, digits=3),
                "total_alt_reads": maybe_num(total_alt, digits=3),
                "expected_alt_fraction": maybe_num(expected_frac, digits=6),
                "blacklist_reason": reason,
                "blacklisted": str(blacklisted),
            })

    return rows


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_uploaded_variation_list(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    values = sorted({
        row["uploaded_variation"]
        for row in rows
        if row.get("uploaded_variation") and row["uploaded_variation"] not in ("", "NA", ".")
    })
    with open(path, "w") as fh:
        for value in values:
            fh.write(f"{value}\n")


def main():
    ap = argparse.ArgumentParser(description="Generate strand-aware blacklist for RNA_EDIT variants from phased VCF, VEP, and RNA BAM.")
    ap.add_argument("--input", action="append", required=True, help="PATIENT=VCF=RNA_BAM")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cohort-support-outfile", default="")
    ap.add_argument("--cohort-blacklist-outfile", default="")
    ap.add_argument("--cohort-uploaded-variation-outfile", default="")
    ap.add_argument("--vep-dir", default="")
    ap.add_argument("--vep-suffix", default="_vep.vep")
    ap.add_argument("--protocol", default="fr-firststrand")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-expected-frac", type=float, default=0.8)
    args = ap.parse_args()

    pairs = []
    for item in args.input:
        parts = item.split("=", 2)
        if len(parts) != 3:
            raise SystemExit(f"--input must be PATIENT=VCF=RNA_BAM, got: {item}")
        patient, vcf_path, bam_path = parts
        pairs.append((patient, vcf_path, bam_path))

    protocol = protocol_alias(args.protocol)
    fieldnames = [
        "SAMPLE",
        "uploaded_variation",
        "location",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene_id",
        "gene_symbol",
        "feature",
        "feature_type",
        "consequence",
        "existing_variation",
        "source_set",
        "filter",
        "protocol",
        "transcript_strand",
        "transcript_strand_set",
        "expected_alt_reads",
        "opposite_alt_reads",
        "total_alt_reads",
        "expected_alt_fraction",
        "blacklist_reason",
        "blacklisted",
    ]

    all_rows = []
    for patient, vcf_path, bam_path in pairs:
        if not os.path.exists(vcf_path):
            print(f"[warn] missing phased VCF for {patient}: {vcf_path}")
            continue
        if not os.path.exists(bam_path):
            print(f"[warn] missing RNA BAM for {patient}: {bam_path}")
            continue
        vep_path = os.path.join(args.vep_dir, f"{patient}{args.vep_suffix}") if args.vep_dir else ""
        if args.vep_dir and not os.path.exists(vep_path):
            print(f"[warn] missing VEP for {patient}: {vep_path}")
            continue

        rows = analyze_patient(
            patient=patient,
            vcf_path=vcf_path,
            bam_path=bam_path,
            vep_path=vep_path,
            protocol=protocol,
            min_mapq=args.min_mapq,
            min_baseq=args.min_baseq,
            min_expected_frac=args.min_expected_frac,
        )
        rows.sort(key=lambda r: (r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
        black_rows = [r for r in rows if r["blacklisted"] == "1"]
        support_path = os.path.join(args.outdir, f"{patient}.strand_support.tsv")
        blacklist_path = os.path.join(args.outdir, f"{patient}.strand_blacklist.tsv")
        uploaded_variation_path = os.path.join(args.outdir, f"{patient}.strand_blacklist.uploaded_variation.txt")
        write_tsv(support_path, rows, fieldnames)
        write_tsv(blacklist_path, black_rows, fieldnames)
        write_uploaded_variation_list(uploaded_variation_path, black_rows)
        print(f"[done] {patient}: support={support_path} blacklist={blacklist_path} uploaded_variation={uploaded_variation_path}")
        all_rows.extend(rows)

    if args.cohort_support_outfile:
        all_rows.sort(key=lambda r: (r["SAMPLE"], r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
        write_tsv(args.cohort_support_outfile, all_rows, fieldnames)
        print(f"[done] cohort support: {args.cohort_support_outfile}")
    if args.cohort_blacklist_outfile:
        black = [r for r in all_rows if r["blacklisted"] == "1"]
        black.sort(key=lambda r: (r["SAMPLE"], r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
        write_tsv(args.cohort_blacklist_outfile, black, fieldnames)
        print(f"[done] cohort blacklist: {args.cohort_blacklist_outfile}")
        if args.cohort_uploaded_variation_outfile:
            write_uploaded_variation_list(args.cohort_uploaded_variation_outfile, black)
            print(f"[done] cohort uploaded_variation blacklist: {args.cohort_uploaded_variation_outfile}")


if __name__ == "__main__":
    main()
