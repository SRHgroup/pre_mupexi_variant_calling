#!/usr/bin/env python3
import argparse
import csv
import os

from annotate_dedup_vep import (
    annotate_rows,
    deduplicate_rows,
    load_vep,
    load_vcf_annotations,
    parse_uploaded_variation,
)


FIELDNAMES = [
    "patient_id",
    "event_type",
    "event_id",
    "mutation_id_vep",
    "uploaded_variation",
    "source_set",
    "known_rnaedit_db",
    "known_db_hit",
    "edit_sig",
    "vcf_filter",
    "Hugo_Symbol",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "primary_consequence",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Location",
    "Allele",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
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


def parse_map(items):
    out = {}
    for item in items:
        if "=" not in item:
            raise SystemExit(f"ERROR: expected PATIENT=/path form, got: {item}")
        patient, path = item.split("=", 1)
        out[patient] = path
    return out


def format_end_position(start_pos, ref):
    if start_pos in (None, "NA"):
        return "NA"
    if not ref or ref == "-":
        return str(start_pos)
    return str(start_pos + max(len(ref), 1) - 1)


def parse_uploaded_variation_extended(uploaded):
    text = str(uploaded).strip()
    if not text or "_" not in text or "/" not in text:
        return None
    try:
        left, allele_change = text.rsplit("_", 1)
        chrom, pos_text = left.rsplit("_", 1)
        ref, alt = allele_change.split("/", 1)
    except ValueError:
        return None
    if not chrom or not pos_text:
        return None
    if "-" in pos_text:
        start_text, end_text = pos_text.split("-", 1)
    else:
        start_text, end_text = pos_text, pos_text
    try:
        start_pos = int(start_text)
        end_pos = int(end_text)
    except Exception:
        return None
    return chrom, pos_text, start_pos, end_pos, ref, alt


def build_mutation_id_vep(uploaded_variation, amino_acids):
    uploaded = str(uploaded_variation).strip()
    if not uploaded or "_" not in uploaded:
        return ""
    if not amino_acids or amino_acids in {"-", "NA"} or "/" not in amino_acids:
        return ""
    aa_normal, aa_mut = str(amino_acids).split("/", 1)
    aa_normal = aa_normal.strip()
    aa_mut = aa_mut.strip()
    if not aa_normal or not aa_mut:
        return ""
    prefix, _alleles = uploaded.rsplit("_", 1)
    return f"{prefix}_{aa_normal}/{aa_mut}"


def build_row(patient, row):
    uploaded_variation = row.get("Uploaded_variation", "NA") or "NA"
    parsed = parse_uploaded_variation_extended(uploaded_variation)
    chrom = "NA"
    pos = "NA"
    end_pos = "NA"
    ref = "NA"
    alt = row.get("Allele", "NA") or "NA"
    if parsed is not None:
        chrom_v, pos_text, pos_v, end_v, ref_v, alt_v = parsed
        chrom = chrom_v
        pos = str(pos_v)
        end_pos = str(end_v)
        ref = ref_v or "NA"
        alt = alt_v or alt
    amino_acids = row.get("Amino_acids", "NA") or "NA"
    mutation_id_vep = build_mutation_id_vep(uploaded_variation, amino_acids)
    return {
        "patient_id": patient,
        "event_type": "SNV",
        "event_id": mutation_id_vep or "NA",
        "mutation_id_vep": mutation_id_vep or "NA",
        "uploaded_variation": uploaded_variation,
        "source_set": row.get("source_set", "NA") or "NA",
        "known_rnaedit_db": row.get("known_rnaedit_db", "NA") or "NA",
        "known_db_hit": row.get("known_db_hit", "0") or "0",
        "edit_sig": row.get("edit_sig", "NA") or "NA",
        "vcf_filter": row.get("vcf_filter", "NA") or "NA",
        "Hugo_Symbol": row.get("gene_symbol", "NA") or "NA",
        "Gene": row.get("Gene", "NA") or "NA",
        "Feature": row.get("Feature", "NA") or "NA",
        "Feature_type": row.get("Feature_type", "NA") or "NA",
        "Consequence": row.get("Consequence", "NA") or "NA",
        "primary_consequence": row.get("primary_consequence", "NA") or "NA",
        "Chromosome": chrom,
        "Start_Position": pos,
        "End_Position": end_pos if end_pos != "NA" else format_end_position(int(pos) if pos != "NA" else None, ref),
        "Reference_Allele": ref,
        "Tumor_Seq_Allele2": alt,
        "Location": row.get("Location", "NA") or "NA",
        "Allele": row.get("Allele", "NA") or "NA",
        "Protein_position": row.get("Protein_position", "NA") or "NA",
        "Amino_acids": amino_acids,
        "Codons": row.get("Codons", "NA") or "NA",
        "Existing_variation": row.get("Existing_variation", "NA") or "NA",
        "gt": row.get("gt", "NA") or "NA",
        "ps": row.get("ps", "NA") or "NA",
        "normal_vaf": row.get("normal_vaf", "NA") or "NA",
        "tumor_vaf": row.get("tumor_vaf", "NA") or "NA",
        "n_alt_count": row.get("n_alt_count", "NA") or "NA",
        "t_alt_count": row.get("t_alt_count", "NA") or "NA",
        "n_ref_count": row.get("n_ref_count", "NA") or "NA",
        "t_ref_count": row.get("t_ref_count", "NA") or "NA",
        "n_depth": row.get("n_depth", "NA") or "NA",
        "t_depth": row.get("t_depth", "NA") or "NA",
        "canonical": row.get("canonical", "0") or "0",
        "has_flags": row.get("has_flags", "0") or "0",
        "is_nmd": row.get("is_nmd", "0") or "0",
        "has_protein": row.get("has_protein", "0") or "0",
        "transcript_row_count": row.get("transcript_row_count", "NA") or "NA",
        "amino_acid_option_count": row.get("amino_acid_option_count", "NA") or "NA",
    }


def load_patient_rows(patient, vep_path, vcf_path, tumor_sample, tumor_labels, normal_labels):
    _meta, _header, rows = load_vep(vep_path)
    by_full, by_alt = load_vcf_annotations(vcf_path, tumor_sample, tumor_labels, normal_labels)
    annotated, _missing_vcf = annotate_rows(rows, by_full, by_alt)
    for row in annotated:
        row["SAMPLE"] = patient
    kept, _filtered, _multi_variant_count = deduplicate_rows(annotated)
    out = []
    for row in kept:
        if row.get("source_set") not in {"SOMATIC", "RNA_EDIT"}:
            continue
        out.append(build_row(patient, row))
    return out


def write_tsv(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in FIELDNAMES})


def main():
    ap = argparse.ArgumentParser(description="Generate a cohort MAF-like summary from MuPeXI VEP files and phased VCFs.")
    ap.add_argument("--vep-input", action="append", required=True, help="PATIENT=/path/to/patient_vep.vep(.gz)")
    ap.add_argument("--vcf-input", action="append", required=True, help="PATIENT=/path/to/phased.vcf.gz")
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--tumor-sample", action="append", default=[], help="PATIENT=TUMOR_SAMPLE_NAME")
    ap.add_argument("--tumor-label", action="append", default=[])
    ap.add_argument("--normal-label", action="append", default=[])
    args = ap.parse_args()

    vep_inputs = parse_map(args.vep_input)
    vcf_inputs = parse_map(args.vcf_input)
    tumor_samples = parse_map(args.tumor_sample)
    tumor_labels = args.tumor_label or ["TUMOR", "DNA_TUMOR", "RNA_TUMOR"]
    normal_labels = args.normal_label or ["DNA_NORMAL"]

    patients = sorted(set(vep_inputs) & set(vcf_inputs))
    if not patients:
        raise SystemExit("ERROR: no overlapping patient VEP/VCF inputs found")

    rows = []
    for patient in patients:
        rows.extend(
            load_patient_rows(
                patient,
                vep_inputs[patient],
                vcf_inputs[patient],
                tumor_samples.get(patient, ""),
                tumor_labels,
                normal_labels,
            )
        )

    rows.sort(key=lambda r: (r["patient_id"], r["event_id"], r["Hugo_Symbol"], r["Protein_position"]))
    write_tsv(args.outfile, rows)
    print(f"[done] wrote {len(rows)} MAF-like rows -> {args.outfile}")


if __name__ == "__main__":
    main()
