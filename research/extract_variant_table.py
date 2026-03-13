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


def parse_info(info):
    out = {}
    if info in (".", ""):
        return out
    for token in info.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
        else:
            out[token] = True
    return out


def parse_format_value(fmt_keys, sample_value, key):
    vals = sample_value.split(":")
    idx = {k: i for i, k in enumerate(fmt_keys)}
    if key not in idx:
        return None
    i = idx[key]
    if i >= len(vals):
        return None
    value = vals[i]
    if value in (".", ""):
        return None
    return value


def parse_numeric(value):
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


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
    if ad:
        parts = ad.split(",")
        wanted = alt_index + 1
        if len(parts) > wanted:
            try:
                return float(parts[wanted])
            except Exception:
                pass
    if fad:
        parts = fad.split(",")
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


def normalize_chrom(chrom):
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


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
    patterns = []
    for label in normal_labels:
        patterns.append(re.compile(rf"(?:_|^){re.escape(label)}(?:\.|$)", re.IGNORECASE))
    for sample in samples:
        is_normal = any(pattern.search(sample) for pattern in patterns)
        if not is_normal:
            return sample
    return None


def parse_vep_location(loc):
    s = str(loc).strip()
    if ":" not in s:
        return None, None
    chrom, rest = s.split(":", 1)
    pos_s = rest.split("-", 1)[0]
    try:
        pos = int(pos_s)
    except Exception:
        return None, None
    return normalize_chrom(chrom), pos


def parse_symbol(extra):
    m = re.search(r"(?:^|;)SYMBOL=([^;]+)", str(extra))
    return m.group(1) if m else ""


def split_amino_acids(value):
    text = str(value).strip()
    if not text or text in (".", "-"):
        return "", ""
    if "/" in text:
        ref, alt = text.split("/", 1)
        return ref, alt
    return text, ""


def add_unique(store, field, value):
    text = str(value).strip()
    if text and text != "." and text != "-":
        store[field].add(text)


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


def fallback_uploaded_variation(chrom, pos, ref, alt):
    return f"{normalize_chrom(chrom)}_{pos}_{ref}/{alt}"


def fallback_location(chrom, pos):
    return f"{normalize_chrom(chrom)}:{pos}"


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
            allele = str(row.get("Allele", "")).strip()
            if chrom is None or pos is None or not allele:
                continue

            key = (chrom, pos, allele)
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
            add_unique(rec, "protein_position", row.get("Protein_position", ""))
            add_unique(rec, "amino_acids", row.get("Amino_acids", ""))
            add_unique(rec, "codons", row.get("Codons", ""))
            add_unique(rec, "existing_variation", row.get("Existing_variation", ""))

            aa_ref, aa_alt = split_amino_acids(row.get("Amino_acids", ""))
            add_unique(rec, "protein_ref", aa_ref)
            add_unique(rec, "protein_alt", aa_alt)

    return out


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
    gt = parse_format_value(fmt_keys, sample_value, "GT")
    ps = parse_format_value(fmt_keys, sample_value, "PS")
    ad = parse_format_value(fmt_keys, sample_value, "AD")
    fad = parse_format_value(fmt_keys, sample_value, "FAD")
    dp = parse_numeric(parse_format_value(fmt_keys, sample_value, "DP"))
    if dp is None:
        dp = infer_dp_from_fields(ad, fad)
    af = parse_alt_specific_numeric(parse_format_value(fmt_keys, sample_value, "AF"), alt_index)
    if af is None:
        af = infer_af_from_ad(ad, alt_index)
    alt_count = infer_alt_count(ad, fad, af, dp, alt_index)
    return {
        "gt": gt or "NA",
        "ps": ps or "NA",
        "dp": maybe_num(dp, digits=3),
        "af": maybe_num(af, digits=6),
        "alt_count": maybe_num(alt_count, digits=3),
    }


def classify_source_set(source_set):
    tokens = [t.strip() for t in str(source_set).split(",") if t.strip()]
    if set(tokens) == {"GERMLINE"}:
        return "GERMLINE"
    if set(tokens) == {"SOMATIC"}:
        return "SOMATIC"
    if set(tokens) == {"RNA_EDIT"}:
        return "RNA_EDIT"
    return str(source_set).strip() or "NA"


def rows_for_patient(patient, vcf_path, vep_map, tumor_sample, tumor_labels, normal_labels):
    rows = []
    with open_text(vcf_path) as fh:
        normal_idx = None
        tumor_idx = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                normal_name = detect_label_sample(samples, "", normal_labels)
                tumor_name = detect_label_sample(samples, tumor_sample, tumor_labels)
                if tumor_name is None:
                    tumor_name = detect_non_normal_sample(samples, normal_labels)
                if normal_name is None or tumor_name is None:
                    raise SystemExit(f"ERROR: could not detect normal/tumor columns in {vcf_path}")
                normal_idx = 9 + samples.index(normal_name)
                tumor_idx = 9 + samples.index(tumor_name)
                continue

            if normal_idx is None or tumor_idx is None:
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < max(normal_idx, tumor_idx) + 1:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")
            flt = fields[6] if fields[6] else "NA"
            info_map = parse_info(fields[7])
            source_set = classify_source_set(info_map.get("SOURCE_SET", ""))
            fmt_keys = fields[8].split(":")

            for alt_index, alt in enumerate(alts):
                normal_metrics = extract_sample_metrics(fmt_keys, fields[normal_idx], alt_index)
                tumor_metrics = extract_sample_metrics(fmt_keys, fields[tumor_idx], alt_index)

                if source_set == "GERMLINE":
                    tumor_metrics["dp"] = "NA"
                    tumor_metrics["af"] = "NA"
                    tumor_metrics["alt_count"] = "NA"

                key = (normalize_chrom(chrom), pos, alt)
                ann = vep_map.get(key, {})
                row = {
                    "SAMPLE": patient,
                    "uploaded_variation": join_values(ann.get("uploaded_variation", {fallback_uploaded_variation(chrom, pos, ref, alt)})),
                    "location": join_values(ann.get("location", {fallback_location(chrom, pos)})),
                    "chrom": normalize_chrom(chrom),
                    "pos": str(pos),
                    "ref": ref,
                    "alt": alt,
                    "protein_ref": join_values(ann.get("protein_ref", set())),
                    "protein_alt": join_values(ann.get("protein_alt", set())),
                    "protein_position": join_values(ann.get("protein_position", set())),
                    "amino_acids": join_values(ann.get("amino_acids", set())),
                    "codons": join_values(ann.get("codons", set())),
                    "gene_id": join_values(ann.get("gene_id", set())),
                    "gene_symbol": join_values(ann.get("gene_symbol", set())),
                    "feature": join_values(ann.get("feature", set())),
                    "feature_type": join_values(ann.get("feature_type", set())),
                    "consequence": join_values(ann.get("consequence", set())),
                    "existing_variation": join_values(ann.get("existing_variation", set())),
                    "source_set": source_set,
                    "filter": flt,
                    "normal_gt": normal_metrics["gt"],
                    "tumor_gt": tumor_metrics["gt"],
                    "tumor_ps": tumor_metrics["ps"],
                    "normal_dp": normal_metrics["dp"],
                    "normal_alt_count": normal_metrics["alt_count"],
                    "normal_vaf": normal_metrics["af"],
                    "tumor_dp": tumor_metrics["dp"],
                    "tumor_alt_count": tumor_metrics["alt_count"],
                    "tumor_vaf": tumor_metrics["af"],
                }
                rows.append(row)
    return rows


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    ap = argparse.ArgumentParser(description="Create a mutation-level research table by joining phased VCF and VEP output.")
    ap.add_argument("--input", action="append", required=True, help="PATIENT=/path/to/phased.vcf.gz")
    ap.add_argument("--outdir", required=True, help="Output directory for per-patient tables")
    ap.add_argument("--cohort-outfile", default="", help="Optional merged cohort TSV path")
    ap.add_argument("--tumor-sample", default="", help="Optional exact tumor sample column name")
    ap.add_argument("--tumor-label", action="append", default=["TUMOR", "DNA_TUMOR", "DNA_TUMOUR"])
    ap.add_argument("--normal-label", action="append", default=["DNA_NORMAL"])
    ap.add_argument("--vep-dir", default="", help="Dir with per-patient VEP files: {patient}_vep.vep")
    ap.add_argument("--vep-suffix", default="_vep.vep", help="Per-patient VEP suffix")
    args = ap.parse_args()

    pairs = []
    for item in args.input:
        if "=" not in item:
            raise SystemExit(f"--input must be PATIENT=VCF, got: {item}")
        patient, vcf = item.split("=", 1)
        pairs.append((patient, vcf))

    fieldnames = [
        "SAMPLE",
        "uploaded_variation",
        "location",
        "chrom",
        "pos",
        "ref",
        "alt",
        "protein_ref",
        "protein_alt",
        "protein_position",
        "amino_acids",
        "codons",
        "gene_id",
        "gene_symbol",
        "feature",
        "feature_type",
        "consequence",
        "existing_variation",
        "source_set",
        "filter",
        "normal_gt",
        "tumor_gt",
        "tumor_ps",
        "normal_dp",
        "normal_alt_count",
        "normal_vaf",
        "tumor_dp",
        "tumor_alt_count",
        "tumor_vaf",
    ]

    all_rows = []
    for patient, vcf_path in pairs:
        if not os.path.exists(vcf_path):
            print(f"[warn] missing phased VCF for {patient}: {vcf_path}")
            continue

        vep_map = {}
        if args.vep_dir:
            vep_path = os.path.join(args.vep_dir, f"{patient}{args.vep_suffix}")
            if os.path.exists(vep_path):
                vep_map = load_vep_map(vep_path)
            else:
                print(f"[warn] missing VEP file for {patient}: {vep_path}")

        rows = rows_for_patient(
            patient=patient,
            vcf_path=vcf_path,
            vep_map=vep_map,
            tumor_sample=args.tumor_sample,
            tumor_labels=args.tumor_label,
            normal_labels=args.normal_label,
        )
        rows.sort(key=lambda r: (r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
        patient_out = os.path.join(args.outdir, f"{patient}.variant_table.tsv")
        write_tsv(patient_out, rows, fieldnames)
        print(f"[done] {patient}: {len(rows)} rows -> {patient_out}")
        all_rows.extend(rows)

    if args.cohort_outfile:
        all_rows.sort(key=lambda r: (r["SAMPLE"], r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
        write_tsv(args.cohort_outfile, all_rows, fieldnames)
        print(f"[done] cohort: {len(all_rows)} rows -> {args.cohort_outfile}")


if __name__ == "__main__":
    main()
