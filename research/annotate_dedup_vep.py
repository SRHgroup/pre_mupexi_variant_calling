#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
from collections import Counter, defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_chrom(chrom):
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def parse_info(info):
    out = {}
    if info in ("", "."):
        return out
    for token in info.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
        elif token:
            out[token] = True
    return out


def parse_format_value(fmt_keys, sample_value, key):
    values = sample_value.split(":")
    idx = {k: i for i, k in enumerate(fmt_keys)}
    if key not in idx:
        return None
    i = idx[key]
    if i >= len(values):
        return None
    value = values[i]
    if value in ("", "."):
        return None
    return value


def parse_numeric(value):
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


def maybe_num(value, digits=6):
    if value is None:
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def infer_af_from_ad(ad, alt_index):
    if not ad:
        return None
    parts = ad.split(",")
    wanted = alt_index + 1
    if len(parts) <= wanted:
        return None
    try:
        ref = float(parts[0])
        alt = float(parts[wanted])
    except Exception:
        return None
    denom = ref + alt
    if denom <= 0:
        return None
    return alt / denom


def infer_alt_count(ad, fad, af, dp, alt_index):
    for value in (ad, fad):
        if not value:
            continue
        parts = value.split(",")
        wanted = alt_index + 1
        if len(parts) > wanted:
            try:
                return float(parts[wanted])
            except Exception:
                pass
    if af is not None and dp is not None:
        try:
            return float(af) * float(dp)
        except Exception:
            return None
    return None


def infer_ref_count(ad, fad):
    for value in (ad, fad):
        if not value:
            continue
        parts = value.split(",")
        if parts:
            try:
                return float(parts[0])
            except Exception:
                pass
    return None


def infer_dp_from_fields(ad, fad):
    for value in (ad, fad):
        if not value:
            continue
        parts = value.split(",")
        try:
            nums = [float(p) for p in parts if p not in ("", ".")]
        except Exception:
            continue
        if nums:
            return sum(nums)
    return None


def parse_alt_specific_numeric(value, alt_index):
    if value is None:
        return None
    parts = str(value).split(",")
    if len(parts) == 1:
        return parse_numeric(parts[0])
    if alt_index < len(parts):
        return parse_numeric(parts[alt_index])
    return None


def extract_sample_metrics(fmt_keys, sample_value, alt_index):
    gt = parse_format_value(fmt_keys, sample_value, "GT") or "NA"
    ps = parse_format_value(fmt_keys, sample_value, "PS") or "NA"
    ad = parse_format_value(fmt_keys, sample_value, "AD")
    fad = parse_format_value(fmt_keys, sample_value, "FAD")
    dp = parse_numeric(parse_format_value(fmt_keys, sample_value, "DP"))
    if dp is None:
        dp = infer_dp_from_fields(ad, fad)
    af = parse_alt_specific_numeric(parse_format_value(fmt_keys, sample_value, "AF"), alt_index)
    if af is None:
        af = infer_af_from_ad(ad, alt_index)
    alt_count = infer_alt_count(ad, fad, af, dp, alt_index)
    ref_count = infer_ref_count(ad, fad)
    return {
        "gt": gt,
        "ps": ps,
        "depth": maybe_num(dp, digits=3),
        "vaf": maybe_num(af, digits=6),
        "alt_count": maybe_num(alt_count, digits=3),
        "ref_count": maybe_num(ref_count, digits=3),
    }


def detect_label_sample(samples, requested, labels):
    if requested:
        for sample in samples:
            if sample == requested:
                return sample
    patterns = []
    for label in labels:
        patterns.append(re.compile(rf"(?:_|^){re.escape(label)}(?:\.|$)", re.IGNORECASE))
        if "TUMOR" in label:
            alt = label.replace("TUMOR", "TUMOUR")
            patterns.append(re.compile(rf"(?:_|^){re.escape(alt)}(?:\.|$)", re.IGNORECASE))
        if "TUMOUR" in label:
            alt = label.replace("TUMOUR", "TUMOR")
            patterns.append(re.compile(rf"(?:_|^){re.escape(alt)}(?:\.|$)", re.IGNORECASE))
    for sample in samples:
        for pattern in patterns:
            if pattern.search(sample):
                return sample
    return None


def detect_non_normal_sample(samples, normal_labels):
    patterns = [re.compile(rf"(?:_|^){re.escape(label)}(?:\.|$)", re.IGNORECASE) for label in normal_labels]
    for sample in samples:
        if not any(pattern.search(sample) for pattern in patterns):
            return sample
    return None


def classify_source_set(value):
    tokens = [t.strip() for t in str(value).split(",") if t.strip()]
    if set(tokens) == {"GERMLINE"}:
        return "GERMLINE"
    if set(tokens) == {"SOMATIC"}:
        return "SOMATIC"
    if set(tokens) == {"RNA_EDIT"}:
        return "RNA_EDIT"
    return str(value).strip() or "NA"


def parse_uploaded_variation(uploaded):
    text = str(uploaded).strip()
    m = re.match(r"(.+?)_(\d+)_([^/]+)/([^/]+)$", text)
    if not m:
        return None
    chrom, pos, ref, alt = m.groups()
    try:
        pos_i = int(pos)
    except Exception:
        return None
    return normalize_chrom(chrom), pos_i, ref, alt


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
    return m.group(1) if m else "NA"


def parse_canonical(extra):
    m = re.search(r"(?:^|;)CANONICAL=([^;]+)", str(extra))
    return "1" if m and m.group(1).upper() in ("YES", "1") else "0"


def has_flags(extra):
    return "1" if re.search(r"(?:^|;)FLAGS=", str(extra)) else "0"


def consequence_rank(value):
    text = str(value)
    if re.search(r"(?:^|,)frameshift_variant(?:,|$)", text):
        return 1, "frameshift_variant"
    if re.search(r"(?:^|,)(?:inframe_insertion|inframe_deletion)(?:,|$)", text):
        if "inframe_insertion" in text:
            return 2, "inframe_insertion"
        return 2, "inframe_deletion"
    if re.search(r"(?:^|,)missense_variant(?:,|$)", text):
        return 3, "missense_variant"
    if re.search(r"(?:^|,)synonymous_variant(?:,|$)", text):
        return 4, "synonymous_variant"
    tokens = [t.strip() for t in text.split(",") if t.strip()]
    return 9, (tokens[0] if tokens else "NA")


def is_nmd(value):
    return "1" if re.search(r"(?:^|,)NMD_transcript_variant(?:,|$)", str(value)) else "0"


def has_protein_fields(row):
    return "1" if str(row.get("Protein_position", "")).strip() not in ("", ".", "-") and str(row.get("Amino_acids", "")).strip() not in ("", ".", "-") else "0"


def load_vep(path):
    meta = []
    header_raw = []
    header = None
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                meta.append(line)
                continue
            if line.startswith("#"):
                header_raw = line.rstrip("\n").split("\t")
                header = [c.lstrip("#") for c in header_raw]
                continue
            if not header:
                continue
            cols = line.rstrip("\n").split("\t")
            row = {header[i]: cols[i] if i < len(cols) else "" for i in range(len(header))}
            row["_raw_cols"] = cols
            rows.append(row)
    return meta, header_raw, rows


def load_vcf_annotations(path, tumor_sample, tumor_labels, normal_labels):
    by_full = {}
    by_alt = {}
    normal_idx = None
    tumor_idx = None
    with open_text(path) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                if samples:
                    normal_name = detect_label_sample(samples, "", normal_labels)
                    tumor_name = detect_label_sample(samples, tumor_sample, tumor_labels)
                    if tumor_name is None:
                        tumor_name = detect_non_normal_sample(samples, normal_labels)
                    if normal_name is not None:
                        normal_idx = 9 + samples.index(normal_name)
                    if tumor_name is not None:
                        tumor_idx = 9 + samples.index(tumor_name)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = normalize_chrom(cols[0])
            try:
                pos = int(cols[1])
            except Exception:
                continue
            ref = cols[3]
            alts = cols[4].split(",")
            flt = cols[6] if cols[6] else "NA"
            info_map = parse_info(cols[7])
            source_set = classify_source_set(info_map.get("SOURCE_SET", ""))
            known_db = str(info_map.get("KNOWN_RNAEDIT_DB", "")).strip() or "NA"
            edit_sig = str(info_map.get("EDIT_SIG", "")).strip() or "NA"
            fmt_keys = cols[8].split(":") if len(cols) > 8 else []
            for alt_index, alt in enumerate(alts):
                normal_metrics = {"gt": "NA", "ps": "NA", "depth": "NA", "vaf": "NA", "alt_count": "NA", "ref_count": "NA"}
                tumor_metrics = {"gt": "NA", "ps": "NA", "depth": "NA", "vaf": "NA", "alt_count": "NA", "ref_count": "NA"}
                if normal_idx is not None and len(cols) > normal_idx:
                    normal_metrics = extract_sample_metrics(fmt_keys, cols[normal_idx], alt_index)
                if tumor_idx is not None and len(cols) > tumor_idx:
                    tumor_metrics = extract_sample_metrics(fmt_keys, cols[tumor_idx], alt_index)
                if source_set == "GERMLINE":
                    tumor_metrics["depth"] = "NA"
                    tumor_metrics["vaf"] = "NA"
                    tumor_metrics["alt_count"] = "NA"
                    tumor_metrics["ref_count"] = "NA"
                rec = {
                    "source_set": source_set,
                    "known_rnaedit_db": known_db,
                    "vcf_filter": flt,
                    "edit_sig": edit_sig,
                    "gt": tumor_metrics["gt"],
                    "ps": tumor_metrics["ps"],
                    "normal_vaf": normal_metrics["vaf"],
                    "tumor_vaf": tumor_metrics["vaf"],
                    "n_alt_count": normal_metrics["alt_count"],
                    "t_alt_count": tumor_metrics["alt_count"],
                    "n_ref_count": normal_metrics["ref_count"],
                    "t_ref_count": tumor_metrics["ref_count"],
                    "n_depth": normal_metrics["depth"],
                    "t_depth": tumor_metrics["depth"],
                }
                by_full[(chrom, pos, ref, alt)] = rec
                by_alt[(chrom, pos, alt)] = rec
    return by_full, by_alt


def annotate_rows(rows, vcf_by_full, vcf_by_alt):
    annotated = []
    missing_vcf = 0
    for row in rows:
        uploaded = str(row.get("Uploaded_variation", "")).strip()
        parsed = parse_uploaded_variation(uploaded)
        if parsed:
            chrom, pos, ref, alt = parsed
        else:
            chrom, pos = parse_vep_location(row.get("Location", ""))
            ref = "NA"
            alt = str(row.get("Allele", "")).strip() or "NA"
        ann = None
        if parsed:
            ann = vcf_by_full.get((chrom, pos, ref, alt))
        if ann is None and chrom is not None and pos is not None and alt not in ("", "NA"):
            ann = vcf_by_alt.get((chrom, pos, alt))
        if ann is None:
            ann = {
                "source_set": "NA",
                "known_rnaedit_db": "NA",
                "vcf_filter": "NA",
                "edit_sig": "NA",
                "gt": "NA",
                "ps": "NA",
                "normal_vaf": "NA",
                "tumor_vaf": "NA",
                "n_alt_count": "NA",
                "t_alt_count": "NA",
                "n_ref_count": "NA",
                "t_ref_count": "NA",
                "n_depth": "NA",
                "t_depth": "NA",
            }
            missing_vcf += 1

        rank, primary = consequence_rank(row.get("Consequence", ""))
        record = dict(row)
        record["SAMPLE"] = "NA"
        record["gene_symbol"] = parse_symbol(row.get("Extra", ""))
        record["canonical"] = parse_canonical(row.get("Extra", ""))
        record["has_flags"] = has_flags(row.get("Extra", ""))
        record["is_nmd"] = is_nmd(row.get("Consequence", ""))
        record["has_protein"] = has_protein_fields(row)
        record["consequence_rank"] = str(rank)
        record["primary_consequence"] = primary
        record["source_set"] = ann["source_set"]
        record["known_rnaedit_db"] = ann["known_rnaedit_db"]
        record["known_db_hit"] = "1" if ann["known_rnaedit_db"] not in ("", ".", "NA") else "0"
        record["vcf_filter"] = ann["vcf_filter"]
        record["edit_sig"] = ann["edit_sig"]
        record["gt"] = ann["gt"]
        record["ps"] = ann["ps"]
        record["normal_vaf"] = ann["normal_vaf"]
        record["tumor_vaf"] = ann["tumor_vaf"]
        record["n_alt_count"] = ann["n_alt_count"]
        record["t_alt_count"] = ann["t_alt_count"]
        record["n_ref_count"] = ann["n_ref_count"]
        record["t_ref_count"] = ann["t_ref_count"]
        record["n_depth"] = ann["n_depth"]
        record["t_depth"] = ann["t_depth"]
        record["dedup_keep"] = "0"
        record["dedup_status"] = "UNSET"
        record["dedup_drop_reason"] = "NA"
        record["transcript_row_count"] = "1"
        record["amino_acid_option_count"] = "1"
        annotated.append(record)
    return annotated, missing_vcf


def deduplicate_rows(rows):
    grouped = defaultdict(list)
    for idx, row in enumerate(rows):
        uploaded = str(row.get("Uploaded_variation", "")).strip()
        key = uploaded if uploaded else f"{row.get('Location', 'NA')}|{row.get('Allele', 'NA')}"
        grouped[key].append((idx, row))

    kept_rows = []
    filtered_rows = []
    multi_variant_count = 0

    for _, group in grouped.items():
        amino_counts = Counter()
        for _, row in group:
            aa = str(row.get("Amino_acids", "")).strip()
            if aa and aa not in (".", "-", "NA"):
                amino_counts[aa] += 1

        aa_options = len(amino_counts) if amino_counts else 1
        if len(group) > 1:
            multi_variant_count += 1

        ranked = []
        for idx, row in group:
            aa = str(row.get("Amino_acids", "")).strip()
            aa_freq = amino_counts.get(aa, 0)
            feature = str(row.get("Feature", "")).strip()
            ranked.append(
                (
                    0 if str(row.get("Feature_type", "")) == "Transcript" else 1,
                    int(row["consequence_rank"]),
                    0 if row["has_protein"] == "1" else 1,
                    0 if row["is_nmd"] == "0" else 1,
                    0 if row["canonical"] == "1" else 1,
                    0 if row["has_flags"] == "0" else 1,
                    -aa_freq,
                    feature,
                    idx,
                    row,
                )
            )
        ranked.sort()
        keep_idx = ranked[0][8]

        for idx, row in group:
            row["transcript_row_count"] = str(len(group))
            row["amino_acid_option_count"] = str(aa_options)
            if idx == keep_idx:
                row["dedup_keep"] = "1"
                row["dedup_status"] = "KEEP_BEST"
                row["dedup_drop_reason"] = "NA"
                kept_rows.append(row)
            else:
                row["dedup_keep"] = "0"
                row["dedup_status"] = "FILTERED_TRANSCRIPT_DUPLICATE"
                row["dedup_drop_reason"] = "LOWER_PRIORITY_TRANSCRIPT"
                filtered_rows.append(row)

    kept_rows.sort(key=lambda r: (r["SAMPLE"], str(r.get("Location", "")), str(r.get("Allele", ""))))
    filtered_rows.sort(key=lambda r: (r["SAMPLE"], str(r.get("Location", "")), str(r.get("Allele", ""))))
    return kept_rows, filtered_rows, multi_variant_count


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def write_stats(path, metrics):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics:
            writer.writerow([key, value])


def build_counts(rows, stage):
    counts = Counter()
    for row in rows:
        counts[
            (
                stage,
                row.get("source_set", "NA"),
                row.get("primary_consequence", "NA"),
                row.get("known_db_hit", "0"),
            )
        ] += 1
    out = []
    for (stage_name, source_set, primary_consequence, known_db_hit), n in sorted(counts.items()):
        out.append(
            {
                "stage": stage_name,
                "source_set": source_set,
                "primary_consequence": primary_consequence,
                "known_db_hit": known_db_hit,
                "count": str(n),
            }
        )
    return out


def patient_fieldnames(vep_header):
    base = [c.lstrip("#") for c in vep_header]
    extra = [
        "SAMPLE",
        "gene_symbol",
        "canonical",
        "has_flags",
        "is_nmd",
        "has_protein",
        "consequence_rank",
        "primary_consequence",
        "source_set",
        "known_rnaedit_db",
        "known_db_hit",
        "vcf_filter",
        "edit_sig",
        "gt",
        "ps",
        "normal_vaf",
        "tumor_vaf",
        "n_alt_count",
        "t_alt_count",
        "n_ref_count",
        "t_ref_count",
        "n_depth",
        "t_depth",
        "dedup_keep",
        "dedup_status",
        "dedup_drop_reason",
        "transcript_row_count",
        "amino_acid_option_count",
    ]
    return base + [c for c in extra if c not in base]


def dedup_fieldnames():
    return [
        "SAMPLE",
        "Uploaded_variation",
        "Location",
        "Allele",
        "Gene",
        "gene_symbol",
        "Feature",
        "Feature_type",
        "Consequence",
        "primary_consequence",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "Existing_variation",
        "source_set",
        "known_rnaedit_db",
        "known_db_hit",
        "vcf_filter",
        "edit_sig",
        "gt",
        "ps",
        "normal_vaf",
        "tumor_vaf",
        "n_alt_count",
        "t_alt_count",
        "n_ref_count",
        "t_ref_count",
        "n_depth",
        "t_depth",
        "canonical",
        "has_flags",
        "is_nmd",
        "has_protein",
        "transcript_row_count",
        "amino_acid_option_count",
    ]


def parse_args():
    ap = argparse.ArgumentParser(description="Merge VCF INFO annotations into VEP rows and deduplicate to one best transcript row per mutation.")
    ap.add_argument("--vep-input", action="append", required=True, help="PATIENT=/path/to/patient_vep.vep")
    ap.add_argument("--vcf-input", action="append", required=True, help="PATIENT=/path/to/phased.vcf.gz")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cohort-dedup-outfile", default="")
    ap.add_argument("--cohort-counts-outfile", default="")
    ap.add_argument("--cohort-stats-outfile", default="")
    ap.add_argument("--tumor-sample", default="")
    ap.add_argument("--tumor-label", action="append", default=["TUMOR", "DNA_TUMOR", "DNA_TUMOUR"])
    ap.add_argument("--normal-label", action="append", default=["DNA_NORMAL"])
    return ap.parse_args()


def parse_input_map(items, label):
    out = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"{label} must be PATIENT=PATH, got: {item}")
        patient, path = item.split("=", 1)
        out[patient] = path
    return out


def row_variant_key(row):
    uploaded = str(row.get("Uploaded_variation", "")).strip()
    if uploaded:
        return uploaded
    return f"{row.get('Location', 'NA')}|{row.get('Allele', 'NA')}"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    vep_inputs = parse_input_map(args.vep_input, "--vep-input")
    vcf_inputs = parse_input_map(args.vcf_input, "--vcf-input")
    patients = sorted(set(vep_inputs) & set(vcf_inputs))
    if not patients:
        raise SystemExit("ERROR: no patients have both VEP and VCF inputs")

    cohort_dedup = []
    cohort_counts = []
    cohort_metric_sums = Counter()

    for patient in patients:
        meta, header_raw, rows = load_vep(vep_inputs[patient])
        vcf_by_full, vcf_by_alt = load_vcf_annotations(
            vcf_inputs[patient],
            tumor_sample=args.tumor_sample,
            tumor_labels=args.tumor_label,
            normal_labels=args.normal_label,
        )
        annotated, missing_vcf = annotate_rows(rows, vcf_by_full, vcf_by_alt)
        for row in annotated:
            row["SAMPLE"] = patient
        kept_rows, filtered_rows, multi_variant_count = deduplicate_rows(annotated)

        unique_variants = {row_variant_key(r) for r in annotated}

        patient_metrics = [
            ("patient", patient),
            ("raw_vep_rows", str(len(annotated))),
            ("raw_unique_variants", str(len(unique_variants))),
            ("variants_with_multiple_transcript_rows", str(multi_variant_count)),
            ("dedup_variants_kept", str(len(kept_rows))),
            ("filtered_transcript_rows", str(len(filtered_rows))),
            ("rows_missing_vcf_annotation", str(missing_vcf)),
        ]

        counts_rows = []
        counts_rows.extend(build_counts(annotated, "raw_rows"))
        counts_rows.extend(build_counts(kept_rows, "dedup_kept_rows"))
        counts_rows.extend(build_counts(filtered_rows, "filtered_transcript_rows"))

        annotated_path = os.path.join(args.outdir, f"{patient}.vep.annotated.tsv")
        dedup_path = os.path.join(args.outdir, f"{patient}.vep.dedup.tsv")
        stats_path = os.path.join(args.outdir, f"{patient}.vep.dedup.stats.tsv")
        counts_path = os.path.join(args.outdir, f"{patient}.vep.dedup.counts.tsv")

        write_tsv(annotated_path, annotated, patient_fieldnames(header_raw))
        write_tsv(dedup_path, kept_rows, dedup_fieldnames())
        write_stats(stats_path, patient_metrics)
        write_tsv(counts_path, counts_rows, ["stage", "source_set", "primary_consequence", "known_db_hit", "count"])

        print(f"[done] {patient}: raw_rows={len(annotated)} dedup_kept={len(kept_rows)} filtered={len(filtered_rows)}")

        cohort_dedup.extend(kept_rows)
        cohort_counts.extend(counts_rows)
        cohort_metric_sums["patients"] += 1
        cohort_metric_sums["raw_vep_rows"] += len(annotated)
        cohort_metric_sums["raw_unique_variants"] += len(unique_variants)
        cohort_metric_sums["variants_with_multiple_transcript_rows"] += multi_variant_count
        cohort_metric_sums["dedup_variants_kept"] += len(kept_rows)
        cohort_metric_sums["filtered_transcript_rows"] += len(filtered_rows)
        cohort_metric_sums["rows_missing_vcf_annotation"] += missing_vcf

    if args.cohort_dedup_outfile:
        cohort_dedup.sort(key=lambda r: (r["SAMPLE"], str(r.get("Location", "")), str(r.get("Allele", ""))))
        write_tsv(args.cohort_dedup_outfile, cohort_dedup, dedup_fieldnames())
        print(f"[done] cohort dedup: {args.cohort_dedup_outfile}")

    if args.cohort_counts_outfile:
        counts_agg = Counter()
        for row in cohort_counts:
            counts_agg[(row["stage"], row["source_set"], row["primary_consequence"], row["known_db_hit"])] += int(row["count"])
        agg_rows = []
        for (stage, source_set, primary_consequence, known_db_hit), count in sorted(counts_agg.items()):
            agg_rows.append(
                {
                    "stage": stage,
                    "source_set": source_set,
                    "primary_consequence": primary_consequence,
                    "known_db_hit": known_db_hit,
                    "count": str(count),
                }
            )
        write_tsv(args.cohort_counts_outfile, agg_rows, ["stage", "source_set", "primary_consequence", "known_db_hit", "count"])
        print(f"[done] cohort counts: {args.cohort_counts_outfile}")

    if args.cohort_stats_outfile:
        write_stats(
            args.cohort_stats_outfile,
            [
                ("patients", str(cohort_metric_sums["patients"])),
                ("raw_vep_rows", str(cohort_metric_sums["raw_vep_rows"])),
                ("raw_unique_variants", str(cohort_metric_sums["raw_unique_variants"])),
                ("variants_with_multiple_transcript_rows", str(cohort_metric_sums["variants_with_multiple_transcript_rows"])),
                ("dedup_variants_kept", str(cohort_metric_sums["dedup_variants_kept"])),
                ("filtered_transcript_rows", str(cohort_metric_sums["filtered_transcript_rows"])),
                ("rows_missing_vcf_annotation", str(cohort_metric_sums["rows_missing_vcf_annotation"])),
            ],
        )
        print(f"[done] cohort stats: {args.cohort_stats_outfile}")


if __name__ == "__main__":
    main()
