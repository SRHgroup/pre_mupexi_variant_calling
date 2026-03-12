#!/usr/bin/env python3
import argparse
import os
import re
from typing import Dict, List, Tuple

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CATEGORY_LABELS = {
    0: "NONE",
    1: "DNA_NORMAL_only",
    2: "DNA_TUMOR_only",
    3: "DNA_NORMAL_DNA_TUMOR",
    4: "RNA_TUMOR_only",
    5: "DNA_NORMAL_RNA_TUMOR",
    6: "DNA_TUMOR_RNA_TUMOR",
    7: "ALL_THREE",
}

STACK_ORDER = [
    "DNA_NORMAL_only",
    "DNA_TUMOR_only",
    "RNA_TUMOR_only",
    "DNA_NORMAL_DNA_TUMOR",
    "DNA_NORMAL_RNA_TUMOR",
    "DNA_TUMOR_RNA_TUMOR",
    "ALL_THREE",
]

STACK_COLORS = {
    "DNA_NORMAL_only": "#4E79A7",
    "DNA_TUMOR_only": "#F28E2B",
    "RNA_TUMOR_only": "#59A14F",
    "DNA_NORMAL_DNA_TUMOR": "#E15759",
    "DNA_NORMAL_RNA_TUMOR": "#76B7B2",
    "DNA_TUMOR_RNA_TUMOR": "#B07AA1",
    "ALL_THREE": "#1F4E79",
}


def sample_base_name(value: str, labels: List[str]) -> str:
    out = value
    for label in labels:
        suffix = f"_{label}"
        if out.endswith(suffix):
            out = out[: -len(suffix)]
    return out


def load_patients(samples_path: str, labels: List[str], only_patient: str) -> List[str]:
    seen = set()
    patients: List[str] = []
    with open(samples_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            sid = re.split(r"[,\t ]+", line)[0]
            if not sid:
                continue
            patient = sample_base_name(sid, labels)
            if not patient:
                continue
            if only_patient and only_patient not in (sid, patient):
                continue
            if patient in seen:
                continue
            seen.add(patient)
            patients.append(patient)
    return patients


def load_md_regions(path: str, cov_col: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", cov_col],
        usecols=[0, 1, 2, 3],
        compression="gzip",
    )


def safe_jaccard(a: pd.Series, b: pd.Series) -> float:
    union = (a | b).sum()
    if union == 0:
        return float("nan")
    return float((a & b).sum() / union)


def compute_patient_overlap(
    patient: str,
    mosdepth_dir: str,
    dna_normal_label: str,
    dna_tumor_label: str,
    rna_tumor_label: str,
    depth_threshold: float,
) -> Tuple[Dict[str, float], pd.DataFrame]:
    dn_path = os.path.join(
        mosdepth_dir,
        f"{patient}_{dna_normal_label}",
        f"{patient}_{dna_normal_label}.md.regions.bed.gz",
    )
    dt_path = os.path.join(
        mosdepth_dir,
        f"{patient}_{dna_tumor_label}",
        f"{patient}_{dna_tumor_label}.md.regions.bed.gz",
    )
    rt_path = os.path.join(
        mosdepth_dir,
        f"{patient}_{rna_tumor_label}",
        f"{patient}_{rna_tumor_label}.md.regions.bed.gz",
    )

    for p in (dn_path, dt_path, rt_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)

    dn = load_md_regions(dn_path, "cov_dn")
    dt = load_md_regions(dt_path, "cov_dt")
    rt = load_md_regions(rt_path, "cov_rt")

    merged = dn.merge(dt, on=["chrom", "start", "end"], how="inner").merge(
        rt, on=["chrom", "start", "end"], how="inner"
    )
    merged["width"] = merged["end"] - merged["start"]

    merged["dn_ok"] = merged["cov_dn"] >= depth_threshold
    merged["dt_ok"] = merged["cov_dt"] >= depth_threshold
    merged["rt_ok"] = merged["cov_rt"] >= depth_threshold
    merged["mask"] = (
        merged["dn_ok"].astype(int)
        + 2 * merged["dt_ok"].astype(int)
        + 4 * merged["rt_ok"].astype(int)
    )
    merged["category"] = merged["mask"].map(CATEGORY_LABELS)

    n_intervals = int(len(merged))
    covered_any = int((merged["mask"] > 0).sum())

    by_cat = (
        merged.groupby("category", as_index=False)
        .agg(n_intervals=("mask", "size"), bp=("width", "sum"))
        .copy()
    )
    by_cat["patient"] = patient
    by_cat["fraction_all"] = by_cat["n_intervals"] / max(n_intervals, 1)
    by_cat["fraction_covered"] = by_cat["n_intervals"] / max(covered_any, 1)

    counts_by_mask = merged["mask"].value_counts().to_dict()
    row: Dict[str, float] = {
        "patient": patient,
        "n_intervals": n_intervals,
        "covered_intervals_any": covered_any,
        "covered_fraction_any": covered_any / max(n_intervals, 1),
        "mean_cov_dna_normal": float(merged["cov_dn"].mean()),
        "mean_cov_dna_tumor": float(merged["cov_dt"].mean()),
        "mean_cov_rna_tumor": float(merged["cov_rt"].mean()),
        "median_cov_dna_normal": float(merged["cov_dn"].median()),
        "median_cov_dna_tumor": float(merged["cov_dt"].median()),
        "median_cov_rna_tumor": float(merged["cov_rt"].median()),
        "jaccard_dn_dt": safe_jaccard(merged["dn_ok"], merged["dt_ok"]),
        "jaccard_dn_rt": safe_jaccard(merged["dn_ok"], merged["rt_ok"]),
        "jaccard_dt_rt": safe_jaccard(merged["dt_ok"], merged["rt_ok"]),
        "pearson_dn_dt": float(merged["cov_dn"].corr(merged["cov_dt"])),
        "pearson_dn_rt": float(merged["cov_dn"].corr(merged["cov_rt"])),
        "pearson_dt_rt": float(merged["cov_dt"].corr(merged["cov_rt"])),
    }
    for mask, label in CATEGORY_LABELS.items():
        row[f"count_{label.lower()}"] = int(counts_by_mask.get(mask, 0))

    return row, by_cat


def plot_stacked_overlap(categories_df: pd.DataFrame, out_png: str, out_svg: str) -> None:
    plot_df = categories_df[categories_df["category"] != "NONE"].copy()
    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.axis("off")
        ax.text(0.5, 0.5, "No covered intervals at selected threshold", ha="center", va="center")
        plt.tight_layout()
        fig.savefig(out_png, dpi=200)
        fig.savefig(out_svg)
        return

    pivot = (
        plot_df.pivot_table(index="patient", columns="category", values="fraction_covered", fill_value=0.0)
        .reindex(columns=STACK_ORDER, fill_value=0.0)
        .copy()
    )
    patients = list(pivot.index)

    fig_h = max(4.0, 0.35 * len(patients) + 2.0)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    left = np.zeros(len(patients))
    y = np.arange(len(patients))
    for cat in STACK_ORDER:
        vals = pivot[cat].values if cat in pivot.columns else np.zeros(len(patients))
        ax.barh(y, vals, left=left, color=STACK_COLORS[cat], edgecolor="white", linewidth=0.3, label=cat)
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels(patients)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of covered intervals (coverage >= threshold)")
    ax.set_title("DNA_NORMAL / DNA_TUMOR / RNA_TUMOR overlap by patient")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=8, frameon=False)
    plt.tight_layout()
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_svg)


def plot_jaccard_heatmap(summary_df: pd.DataFrame, out_png: str, out_svg: str) -> None:
    cols = ["jaccard_dn_dt", "jaccard_dn_rt", "jaccard_dt_rt"]
    labels = ["DN vs DT", "DN vs RT", "DT vs RT"]
    heat = summary_df[cols].to_numpy(dtype=float)
    patients = summary_df["patient"].tolist()

    fig_h = max(4.0, 0.35 * len(patients) + 1.8)
    fig, ax = plt.subplots(figsize=(8, fig_h))
    im = ax.imshow(heat, aspect="auto", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticks(np.arange(len(patients)))
    ax.set_yticklabels(patients)
    ax.set_title("Pairwise Jaccard overlap (coverage >= threshold)")
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Jaccard index")
    plt.tight_layout()
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_svg)


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot overlap between DNA normal, DNA tumor, RNA tumor using mosdepth md.regions.bed.gz")
    ap.add_argument("--samples", required=True, help="SAMPLES.tsv path")
    ap.add_argument("--mosdepth-dir", required=True, help="Directory containing mosdepth patient folders")
    ap.add_argument("--dna-normal-label", default="DNA_NORMAL")
    ap.add_argument("--dna-tumor-label", default="DNA_TUMOR")
    ap.add_argument("--rna-tumor-label", default="RNA_TUMOR")
    ap.add_argument("--depth-threshold", type=float, default=10.0, help="Coverage threshold for overlap")
    ap.add_argument("--patient", default="", help="Optional single patient")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    labels = [
        args.dna_normal_label,
        args.dna_tumor_label,
        args.rna_tumor_label,
        "DNA_NORMAL",
        "DNA_TUMOR",
        "DNA_TUMOUR",
        "RNA_TUMOR",
        "RNA_TUMOUR",
    ]
    patients = load_patients(args.samples, labels, args.patient)
    if not patients:
        raise SystemExit("ERROR: no patients discovered from samples file")

    summary_rows: List[Dict[str, float]] = []
    category_rows: List[pd.DataFrame] = []
    missing_rows: List[Dict[str, str]] = []

    for patient in patients:
        try:
            row, by_cat = compute_patient_overlap(
                patient=patient,
                mosdepth_dir=args.mosdepth_dir,
                dna_normal_label=args.dna_normal_label,
                dna_tumor_label=args.dna_tumor_label,
                rna_tumor_label=args.rna_tumor_label,
                depth_threshold=args.depth_threshold,
            )
            summary_rows.append(row)
            category_rows.append(by_cat)
            print(f"[ok] {patient}")
        except FileNotFoundError as e:
            missing_rows.append({"patient": patient, "missing_file": str(e)})
            print(f"[skip] {patient}: missing {e}")

    if not summary_rows:
        raise SystemExit("ERROR: no patients processed successfully (all missing mosdepth inputs?)")

    summary_df = pd.DataFrame(summary_rows).sort_values("patient").reset_index(drop=True)
    categories_df = pd.concat(category_rows, ignore_index=True)

    summary_tsv = os.path.join(args.outdir, "mosdepth_overlap_summary.tsv")
    categories_tsv = os.path.join(args.outdir, "mosdepth_overlap_categories.tsv")
    missing_tsv = os.path.join(args.outdir, "mosdepth_overlap_missing_inputs.tsv")
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    categories_df.to_csv(categories_tsv, sep="\t", index=False)
    pd.DataFrame(missing_rows).to_csv(missing_tsv, sep="\t", index=False)

    stacked_png = os.path.join(args.outdir, "mosdepth_overlap_stacked.png")
    stacked_svg = os.path.join(args.outdir, "mosdepth_overlap_stacked.svg")
    jaccard_png = os.path.join(args.outdir, "mosdepth_overlap_jaccard_heatmap.png")
    jaccard_svg = os.path.join(args.outdir, "mosdepth_overlap_jaccard_heatmap.svg")

    plot_stacked_overlap(categories_df, stacked_png, stacked_svg)
    plot_jaccard_heatmap(summary_df, jaccard_png, jaccard_svg)

    print(f"[done] summary: {summary_tsv}")
    print(f"[done] categories: {categories_tsv}")
    print(f"[done] missing: {missing_tsv}")
    print(f"[done] plots: {stacked_png}, {stacked_svg}, {jaccard_png}, {jaccard_svg}")


if __name__ == "__main__":
    main()
