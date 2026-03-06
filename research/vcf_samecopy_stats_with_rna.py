#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import sys
from bisect import bisect_left
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

ALLOWED_PROTEIN_CODING_CONSEQUENCES = {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "frameshift_variant",
}


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def file_nonempty(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def count_vcf_records(path: str) -> int:
    n = 0
    with open_text(path) as fh:
        for line in fh:
            if line and not line.startswith("#"):
                n += 1
    return n


def normalize_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def parse_vep_location(loc: str) -> Tuple[Optional[str], Optional[int]]:
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


def load_vep_allowed_variant_keys(path: str) -> set:
    allowed = set()
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            chrom, pos = parse_vep_location(cols[1])
            if chrom is None or pos is None:
                continue
            allele = cols[2].strip()
            consequence = cols[6].strip()
            if not allele or not consequence:
                continue
            # VEP can provide multiple terms separated by comma and/or '&'
            terms = [t for t in re.split(r"[,&]", consequence) if t]
            if not any(t in ALLOWED_PROTEIN_CODING_CONSEQUENCES for t in terms):
                continue
            # If BIOTYPE is present in EXTRA, require protein_coding.
            extra = cols[13] if len(cols) > 13 else ""
            m = re.search(r"(?:^|;)BIOTYPE=([^;]+)", extra)
            if m and m.group(1) != "protein_coding":
                continue
            allowed.add((chrom, pos, allele))
    return allowed


def variant_in_allowed_set(chrom: str, pos: int, alt: str, allowed_keys: set) -> bool:
    c = normalize_chrom(chrom)
    for a in str(alt).split(","):
        aa = a.strip()
        if aa and (c, pos, aa) in allowed_keys:
            return True
    return False


def count_vcf_records_with_allowed(path: str, allowed_keys: set) -> int:
    n = 0
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 5:
                continue
            chrom = toks[0]
            pos = int(toks[1])
            alt = toks[4]
            if variant_in_allowed_set(chrom, pos, alt, allowed_keys):
                n += 1
    return n


def parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if info == "." or not info:
        return out
    for tok in info.split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
        else:
            out[tok] = "1"
    return out


def parse_gt(gt: str) -> Tuple[str, Optional[List[str]]]:
    if not gt or gt == ".":
        return "", None
    if "|" in gt:
        sep = "|"
    elif "/" in gt:
        sep = "/"
    else:
        return "", None
    parts = gt.split(sep)
    if len(parts) != 2:
        return sep, None
    return sep, parts


def get_fmt_index(fmt: str, key: str) -> int:
    parts = fmt.split(":")
    for i, p in enumerate(parts):
        if p == key:
            return i
    return -1


def get_sample_field(sample: str, idx: int) -> str:
    if idx < 0:
        return ""
    parts = sample.split(":")
    if idx >= len(parts):
        return ""
    v = parts[idx]
    return "" if v == "." else v


def is_homo_alt(alleles: List[str]) -> bool:
    if "." in alleles:
        return False
    return alleles[0] == alleles[1] and alleles[0] != "0"


def has_alt(alleles: List[str]) -> bool:
    return any(a not in ("0", ".") for a in alleles)


def alt_sides_if_phased(sep: str, alleles: List[str]) -> set:
    if sep != "|" or "." in alleles:
        return set()
    out = set()
    if alleles[0] not in ("0", "."):
        out.add(0)
    if alleles[1] not in ("0", "."):
        out.add(1)
    return out


def range_slice(sorted_positions: List[int], pos: int, window: int) -> Tuple[int, int]:
    lo = bisect_left(sorted_positions, pos - window)
    hi = bisect_left(sorted_positions, pos + window + 1)
    return lo, hi


def choose_tumor_col(header_samples: List[str], preferred_tumor_label: str) -> int:
    if not header_samples:
        return 9
    # Exact preferred first
    for i, s in enumerate(header_samples):
        if s == preferred_tumor_label:
            return 9 + i
    # Common fallback names
    for i, s in enumerate(header_samples):
        u = s.upper()
        if u == "TUMOR" or u == "TUMOUR":
            return 9 + i
    # Name contains tumor
    for i, s in enumerate(header_samples):
        u = s.upper()
        if "TUMOR" in u or "TUMOUR" in u:
            return 9 + i
    # Two-sample fallback: non-normal
    if len(header_samples) == 2:
        a, b = header_samples[0], header_samples[1]
        au, bu = a.upper(), b.upper()
        a_norm = "NORMAL" in au or au.endswith("_N")
        b_norm = "NORMAL" in bu or bu.endswith("_N")
        if a_norm and not b_norm:
            return 10
        if b_norm and not a_norm:
            return 9
    return 9


def classify_source(info_map: Dict[str, str]) -> str:
    s = str(info_map.get("SOURCE_SET", ""))
    if "GERMLINE" in s:
        return "GERMLINE"
    if "SOMATIC" in s:
        return "SOMATIC"
    if "RNA_EDIT" in s:
        return "RNA_EDIT"
    return "OTHER"


def sample_base_name(value: str, labels: Iterable[str]) -> str:
    out = value
    for lab in labels:
        if not lab:
            continue
        if out.endswith("_" + lab):
            out = out[: -(len(lab) + 1)]
    return out


def discover_patients(samples_path: str, labels: Iterable[str], selected: str) -> List[str]:
    out: List[str] = []
    seen = set()
    with open(samples_path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            sid = re.split(r"[,\t ]+", s)[0]
            pat = sample_base_name(sid, labels)
            if not pat:
                continue
            if selected and selected not in (sid, pat):
                continue
            if pat in seen:
                continue
            seen.add(pat)
            out.append(pat)
    return out


@dataclass(frozen=True)
class VarRec:
    pos: int
    phased_gt: bool
    unphased_gt: bool
    ps: str
    alt_sides: frozenset
    homo_alt: bool
    het_alt: bool


def is_samecopy(a: VarRec, b: VarRec) -> Tuple[bool, bool]:
    if a.homo_alt or b.homo_alt:
        return True, False
    if not (a.phased_gt and b.phased_gt):
        return False, False
    if not a.ps or a.ps == "." or not b.ps or b.ps == ".":
        return False, False
    if a.ps != b.ps:
        return False, False
    if a.alt_sides and b.alt_sides and (set(a.alt_sides) & set(b.alt_sides)):
        return True, True
    return False, False


def first_existing(paths: List[str]) -> str:
    for p in paths:
        if file_nonempty(p):
            return p
    return ""


def resolve_src_suffix(template: str, patient: str) -> str:
    out = template.replace("{patient}", patient)
    return out.rstrip("}")


def parse_merged(merged_vcf: str, tumor_col: int, allowed_somatic_rna_keys: Optional[set] = None):
    by_source: Dict[str, Dict[str, List[VarRec]]] = {
        "GERMLINE": {},
        "SOMATIC": {},
        "RNA_EDIT": {},
        "RNA_EDIT_KNOWN": {},
    }
    counts = {
        "merged_total": 0,
        "merged_germline_total": 0,
        "merged_somatic_total": 0,
        "merged_rna_edit_total": 0,
        "merged_rna_edit_adar_total": 0,
        "merged_rna_edit_apobec3_total": 0,
        "merged_rna_edit_known_db_total": 0,
        "merged_other_total": 0,
        "merged_germline_gt_phased": 0,
        "merged_germline_gt_unphased": 0,
        "merged_somatic_gt_phased": 0,
        "merged_somatic_gt_unphased": 0,
        "merged_rna_edit_gt_phased": 0,
        "merged_rna_edit_gt_unphased": 0,
        "merged_germline_hom_alt": 0,
        "merged_germline_het_alt": 0,
        "merged_somatic_hom_alt": 0,
        "merged_somatic_het_alt": 0,
        "merged_rna_edit_hom_alt": 0,
        "merged_rna_edit_het_alt": 0,
    }

    with open_text(merged_vcf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            counts["merged_total"] += 1
            toks = line.rstrip("\n").split("\t")
            if len(toks) <= max(8, tumor_col):
                counts["merged_other_total"] += 1
                continue
            chrom = toks[0]
            pos = int(toks[1])
            info_map = parse_info(toks[7])
            source = classify_source(info_map)
            if source not in by_source:
                counts["merged_other_total"] += 1
                continue
            if source in ("SOMATIC", "RNA_EDIT") and allowed_somatic_rna_keys is not None:
                alt = toks[4]
                if not variant_in_allowed_set(chrom, pos, alt, allowed_somatic_rna_keys):
                    continue

            fmt = toks[8]
            sample = toks[tumor_col]
            gt_i = get_fmt_index(fmt, "GT")
            ps_i = get_fmt_index(fmt, "PS")
            gt = get_sample_field(sample, gt_i)
            ps = get_sample_field(sample, ps_i)
            sep, alleles = parse_gt(gt)
            phased_gt = sep == "|"
            unphased_gt = sep == "/"
            homo_alt = False
            het_alt = False
            sides = frozenset()
            if alleles is not None and has_alt(alleles):
                homo_alt = is_homo_alt(alleles)
                het_alt = not homo_alt
                if phased_gt:
                    sides = frozenset(alt_sides_if_phased(sep, alleles))

            rec = VarRec(
                pos=pos,
                phased_gt=phased_gt,
                unphased_gt=unphased_gt,
                ps=ps,
                alt_sides=sides,
                homo_alt=homo_alt,
                het_alt=het_alt,
            )
            by_source[source].setdefault(chrom, []).append(rec)

            low = source.lower()
            counts[f"merged_{low}_total"] += 1
            if source == "RNA_EDIT":
                edit_sig = str(info_map.get("EDIT_SIG", "")).upper()
                if "ADAR" in edit_sig:
                    counts["merged_rna_edit_adar_total"] += 1
                if "APOBEC3" in edit_sig:
                    counts["merged_rna_edit_apobec3_total"] += 1
                known_db = str(info_map.get("KNOWN_RNAEDIT_DB", ""))
                if known_db not in ("", "."):
                    counts["merged_rna_edit_known_db_total"] += 1
                    by_source["RNA_EDIT_KNOWN"].setdefault(chrom, []).append(rec)
            if phased_gt:
                counts[f"merged_{low}_gt_phased"] += 1
            if unphased_gt:
                counts[f"merged_{low}_gt_unphased"] += 1
            if homo_alt:
                counts[f"merged_{low}_hom_alt"] += 1
            elif het_alt:
                counts[f"merged_{low}_het_alt"] += 1

    for source in by_source:
        for chrom in by_source[source]:
            by_source[source][chrom].sort(key=lambda r: r.pos)

    by_pos: Dict[str, Dict[str, List[int]]] = {}
    for source in by_source:
        by_pos[source] = {chrom: [r.pos for r in rows] for chrom, rows in by_source[source].items()}

    return by_source, by_pos, counts


def adjacency_counts(src_a: Dict[str, List[VarRec]], pos_b: Dict[str, List[int]], src_b: Dict[str, List[VarRec]], window: int):
    hom_adj_any = 0
    het_samecopy_any = 0
    het_samecopy_phasedcis = 0
    for chrom, rows_a in src_a.items():
        p_b = pos_b.get(chrom)
        rows_b = src_b.get(chrom)
        if not p_b or not rows_b:
            continue
        for a in rows_a:
            lo, hi = range_slice(p_b, a.pos, window)
            if lo == hi:
                continue
            candidates = rows_b[lo:hi]
            if a.homo_alt:
                hom_adj_any += 1
                continue
            if not a.het_alt:
                continue
            found_any = False
            found_cis = False
            for b in candidates:
                same, cis = is_samecopy(a, b)
                if same:
                    found_any = True
                    if cis and b.het_alt:
                        found_cis = True
                    break
            if found_any:
                het_samecopy_any += 1
            if found_cis:
                het_samecopy_phasedcis += 1
    return hom_adj_any, het_samecopy_any, het_samecopy_phasedcis


def main():
    ap = argparse.ArgumentParser(description="Cohort same-copy stats from merged phased VCF including RNA_EDIT.")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--vcfdir", required=True)
    ap.add_argument("--out-normal-label", default="DNA_NORMAL")
    ap.add_argument("--out-dna-label", default="DNA_TUMOR")
    ap.add_argument("--out-rna-label", default="RNA_TUMOR")
    ap.add_argument("--gdna3-ext", required=True)
    ap.add_argument("--source-dna-ext", required=True)
    ap.add_argument("--source-rna-ext", required=True)
    ap.add_argument("--merged-ext", required=True)
    ap.add_argument("--window", type=int, default=33)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--patient", default="")
    ap.add_argument("--tumor-col-label", default="TUMOR")
    ap.add_argument("--vep-dir", default="", help="Dir with per-patient VEP files: {patient}_vep.vep")
    ap.add_argument("--vep-suffix", default="_vep.vep", help="Per-patient VEP suffix")
    args = ap.parse_args()

    labels = [
        "DNA_NORMAL",
        "DNA_TUMOR",
        "DNA_TUMOUR",
        "RNA_TUMOR",
        "RNA_TUMOUR",
        args.out_normal_label,
        args.out_dna_label,
        args.out_rna_label,
    ]
    patients = discover_patients(args.samples, labels, args.patient)
    if not patients:
        print("ERROR: no patients discovered", file=sys.stderr)
        sys.exit(1)

    cols = [
        "patient",
        "germline_selected_total",
        "somatic_filtered_total",
        "rna_filtered_total",
        "merged_total",
        "merged_germline_total",
        "merged_somatic_total",
        "merged_rna_edit_total",
        "merged_rna_edit_adar_total",
        "merged_rna_edit_apobec3_total",
        "merged_rna_edit_known_db_total",
        "merged_other_total",
        "merged_germline_gt_phased",
        "merged_germline_gt_unphased",
        "merged_somatic_gt_phased",
        "merged_somatic_gt_unphased",
        "merged_rna_edit_gt_phased",
        "merged_rna_edit_gt_unphased",
        "merged_germline_hom_alt",
        "merged_germline_het_alt",
        "merged_somatic_hom_alt",
        "merged_somatic_het_alt",
        "merged_rna_edit_hom_alt",
        "merged_rna_edit_het_alt",
        "germline_hom_adj_somatic_any",
        "germline_het_samecopy_with_somatic",
        "germline_het_samecopy_with_somatic_phasedcis",
        "somatic_hom_adj_germline_any",
        "somatic_het_samecopy_with_germline",
        "somatic_het_samecopy_with_germline_phasedcis",
        "germline_hom_adj_rna_any",
        "germline_het_samecopy_with_rna",
        "germline_het_samecopy_with_rna_phasedcis",
        "rna_hom_adj_germline_any",
        "rna_het_samecopy_with_germline",
        "rna_het_samecopy_with_germline_phasedcis",
        "somatic_hom_adj_rna_any",
        "somatic_het_samecopy_with_rna",
        "somatic_het_samecopy_with_rna_phasedcis",
        "rna_hom_adj_somatic_any",
        "rna_het_samecopy_with_somatic",
        "rna_het_samecopy_with_somatic_phasedcis",
    ]

    os.makedirs(os.path.dirname(args.outfile) or ".", exist_ok=True)
    written_rows = 0
    with open(args.outfile, "w") as out:
        out.write("\t".join(cols) + "\n")
        for i, patient in enumerate(patients, 1):
            print(f"[{i}/{len(patients)}] {patient}", file=sys.stderr)
            vep_path = ""
            allowed_keys = None
            if args.vep_dir:
                vep_path = os.path.join(args.vep_dir, f"{patient}{args.vep_suffix}")
                if not file_nonempty(vep_path):
                    print(f"[skip] {patient}: missing VEP file for coding filter: {vep_path}", file=sys.stderr)
                    continue
                allowed_keys = load_vep_allowed_variant_keys(vep_path)
                if not allowed_keys:
                    print(f"[skip] {patient}: no allowed coding SOMATIC/RNA variants in VEP: {vep_path}", file=sys.stderr)
                    continue
            germ = os.path.join(
                args.vcfdir,
                f"{patient}_{args.out_normal_label}",
                f"{patient}_{args.gdna3_ext}",
            )
            if not file_nonempty(germ) and file_nonempty(germ + ".gz"):
                germ = germ + ".gz"
            dna_src_suffix = resolve_src_suffix(args.source_dna_ext, patient)
            rna_src_suffix = resolve_src_suffix(args.source_rna_ext, patient)
            dna_src = first_existing(
                [
                    os.path.join(args.vcfdir, f"{patient}_{args.out_dna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{dna_src_suffix}"),
                    os.path.join(args.vcfdir, f"{patient}_{args.out_dna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{dna_src_suffix}".removesuffix(".gz")),
                    os.path.join(args.vcfdir, f"{patient}_{args.out_dna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{args.out_dna_label}_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf.gz"),
                    os.path.join(args.vcfdir, f"{patient}_{args.out_dna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{args.out_dna_label}_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf"),
                    os.path.join(args.vcfdir, f"{patient}_DNA_TUMOR_vs_{patient}_{args.out_normal_label}", f"{patient}_DNA_TUMOR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf.gz"),
                    os.path.join(args.vcfdir, f"{patient}_DNA_TUMOR_vs_{patient}_{args.out_normal_label}", f"{patient}_DNA_TUMOR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf"),
                    os.path.join(args.vcfdir, f"{patient}_DNA_TUMOUR_vs_{patient}_{args.out_normal_label}", f"{patient}_DNA_TUMOUR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf.gz"),
                    os.path.join(args.vcfdir, f"{patient}_DNA_TUMOUR_vs_{patient}_{args.out_normal_label}", f"{patient}_DNA_TUMOUR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf"),
                ]
            )
            rna_src = first_existing(
                [
                    os.path.join(args.vcfdir, f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{rna_src_suffix}"),
                    os.path.join(args.vcfdir, f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{rna_src_suffix}".removesuffix(".gz")),
                    os.path.join(args.vcfdir, f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf.gz"),
                    os.path.join(args.vcfdir, f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}", f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf"),
                    os.path.join(args.vcfdir, f"{patient}_RNA_TUMOR_vs_{patient}_{args.out_normal_label}", f"{patient}_RNA_TUMOR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf.gz"),
                    os.path.join(args.vcfdir, f"{patient}_RNA_TUMOR_vs_{patient}_{args.out_normal_label}", f"{patient}_RNA_TUMOR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf"),
                    os.path.join(args.vcfdir, f"{patient}_RNA_TUMOUR_vs_{patient}_{args.out_normal_label}", f"{patient}_RNA_TUMOUR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf.gz"),
                    os.path.join(args.vcfdir, f"{patient}_RNA_TUMOUR_vs_{patient}_{args.out_normal_label}", f"{patient}_RNA_TUMOUR_vs_{patient}_{args.out_normal_label}.mutect2.filtered.vcf"),
                ]
            )
            merged = os.path.join(
                args.vcfdir,
                f"{patient}_{args.out_rna_label}_vs_{patient}_{args.out_normal_label}",
                f"{patient}_{args.merged_ext}",
            )
            if not file_nonempty(merged) and file_nonempty(merged + ".gz"):
                merged = merged + ".gz"

            miss = [p for p in [germ, dna_src, rna_src, merged] if not file_nonempty(p)]
            if miss:
                print(f"[skip] {patient}: missing input(s): {', '.join(miss)}", file=sys.stderr)
                continue

            germ_total = count_vcf_records(germ)
            if allowed_keys is not None:
                som_total = count_vcf_records_with_allowed(dna_src, allowed_keys)
                rna_total = count_vcf_records_with_allowed(rna_src, allowed_keys)
            else:
                som_total = count_vcf_records(dna_src)
                rna_total = count_vcf_records(rna_src)

            tumor_col = 9
            with open_text(merged) as fh:
                for line in fh:
                    if line.startswith("#CHROM"):
                        toks = line.rstrip("\n").split("\t")
                        tumor_col = choose_tumor_col(toks[9:], args.tumor_col_label)
                        break

            by_source, by_pos, mc = parse_merged(merged, tumor_col, allowed_keys)

            gs = adjacency_counts(by_source["GERMLINE"], by_pos["SOMATIC"], by_source["SOMATIC"], args.window)
            sg = adjacency_counts(by_source["SOMATIC"], by_pos["GERMLINE"], by_source["GERMLINE"], args.window)
            # For RNA-involving samecopy/phasing metrics, use only KNOWN_RNAEDIT_DB RNA variants.
            gr = adjacency_counts(by_source["GERMLINE"], by_pos["RNA_EDIT_KNOWN"], by_source["RNA_EDIT_KNOWN"], args.window)
            rg = adjacency_counts(by_source["RNA_EDIT_KNOWN"], by_pos["GERMLINE"], by_source["GERMLINE"], args.window)
            sr = adjacency_counts(by_source["SOMATIC"], by_pos["RNA_EDIT_KNOWN"], by_source["RNA_EDIT_KNOWN"], args.window)
            rs = adjacency_counts(by_source["RNA_EDIT_KNOWN"], by_pos["SOMATIC"], by_source["SOMATIC"], args.window)

            row = {
                "patient": patient,
                "germline_selected_total": germ_total,
                "somatic_filtered_total": som_total,
                "rna_filtered_total": rna_total,
                **mc,
                "germline_hom_adj_somatic_any": gs[0],
                "germline_het_samecopy_with_somatic": gs[1],
                "germline_het_samecopy_with_somatic_phasedcis": gs[2],
                "somatic_hom_adj_germline_any": sg[0],
                "somatic_het_samecopy_with_germline": sg[1],
                "somatic_het_samecopy_with_germline_phasedcis": sg[2],
                "germline_hom_adj_rna_any": gr[0],
                "germline_het_samecopy_with_rna": gr[1],
                "germline_het_samecopy_with_rna_phasedcis": gr[2],
                "rna_hom_adj_germline_any": rg[0],
                "rna_het_samecopy_with_germline": rg[1],
                "rna_het_samecopy_with_germline_phasedcis": rg[2],
                "somatic_hom_adj_rna_any": sr[0],
                "somatic_het_samecopy_with_rna": sr[1],
                "somatic_het_samecopy_with_rna_phasedcis": sr[2],
                "rna_hom_adj_somatic_any": rs[0],
                "rna_het_samecopy_with_somatic": rs[1],
                "rna_het_samecopy_with_somatic_phasedcis": rs[2],
            }
            out.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
            written_rows += 1

    if written_rows == 0:
        print(f"ERROR: no patient rows written; check skip messages above. Output header only: {args.outfile}", file=sys.stderr)
        sys.exit(2)
    print(f"[done] wrote {written_rows} patient rows -> {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
