#!/usr/bin/env python3
import argparse
import glob
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


def normalize_chrom(chrom: str) -> str:
    s = str(chrom)
    if s.startswith("chr"):
        s = s[3:]
    return s


def is_canonical_chrom(chrom: str) -> bool:
    s = normalize_chrom(chrom).upper()
    if s in {"X", "Y", "M", "MT"}:
        return True
    return s.isdigit() and 1 <= int(s) <= 22


def chrom_sort_key(chrom: str) -> Tuple[int, int, str]:
    s = normalize_chrom(chrom)
    if s.isdigit():
        return (0, int(s), str(chrom))
    if s.upper() == "X":
        return (1, 23, str(chrom))
    if s.upper() == "Y":
        return (1, 24, str(chrom))
    if s.upper() in ("M", "MT"):
        return (2, 25, str(chrom))
    return (3, 999, str(chrom))


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text))


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


def clear_chromosome_plot_outputs(outdir: str) -> int:
    patterns = [
        "mosdepth_overlap_*_cohort_profile.png",
        "mosdepth_overlap_*_cohort_profile.svg",
        "mosdepth_overlap_*_all_three_overlap_heatmap.png",
        "mosdepth_overlap_*_all_three_overlap_heatmap.svg",
        "mosdepth_overlap_*_discordant_overlap_heatmap.png",
        "mosdepth_overlap_*_discordant_overlap_heatmap.svg",
    ]
    removed = 0
    for pattern in patterns:
        for path in glob.glob(os.path.join(outdir, pattern)):
            try:
                os.remove(path)
                removed += 1
            except OSError:
                pass
    return removed


def compute_patient_overlap(
    patient: str,
    mosdepth_dir: str,
    dna_normal_label: str,
    dna_tumor_label: str,
    rna_tumor_label: str,
    depth_threshold: float,
    region_bin_size: int,
) -> Tuple[Dict[str, float], pd.DataFrame, pd.DataFrame]:
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
    merged = merged[merged["chrom"].map(is_canonical_chrom)].copy()
    if merged.empty:
        raise ValueError(f"{patient}: no canonical chromosomes found in merged mosdepth intervals")
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

    merged["bin_start"] = (merged["start"] // region_bin_size) * region_bin_size
    by_bin = (
        merged.groupby(["chrom", "bin_start"], as_index=False)
        .agg(
            n_intervals=("mask", "size"),
            frac_any=("mask", lambda s: float((s > 0).mean())),
            frac_all_three=("mask", lambda s: float((s == 7).mean())),
            frac_discordant=("mask", lambda s: float(((s > 0) & (s < 7)).mean())),
            mean_cov_dna_normal=("cov_dn", "mean"),
            mean_cov_dna_tumor=("cov_dt", "mean"),
            mean_cov_rna_tumor=("cov_rt", "mean"),
        )
        .copy()
    )
    by_bin["patient"] = patient

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

    return row, by_cat, by_bin


def plot_stacked_overlap(categories_df: pd.DataFrame, out_png: str, out_svg: str) -> None:
    plot_df = categories_df[categories_df["category"] != "NONE"].copy()
    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.axis("off")
        ax.text(0.5, 0.5, "No covered intervals at selected threshold", ha="center", va="center")
        plt.tight_layout()
        fig.savefig(out_png, dpi=200)
        fig.savefig(out_svg)
        plt.close(fig)
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
    plt.close(fig)


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
    plt.close(fig)


def plot_chromosome_profiles(bin_df: pd.DataFrame, outdir: str, region_bin_size: int) -> List[str]:
    produced: List[str] = []
    chroms = sorted(bin_df["chrom"].dropna().unique(), key=chrom_sort_key)
    for chrom in chroms:
        sub = bin_df[bin_df["chrom"] == chrom].copy()
        if sub.empty:
            continue
        agg = (
            sub.groupby("bin_start", as_index=False)
            .agg(
                frac_any=("frac_any", "mean"),
                frac_all_three=("frac_all_three", "mean"),
                frac_discordant=("frac_discordant", "mean"),
            )
            .sort_values("bin_start")
        )
        if agg.empty:
            continue

        x = (agg["bin_start"].to_numpy(dtype=float) + (region_bin_size / 2.0)) / 1e6
        fig, ax = plt.subplots(figsize=(14, 4))
        ax.plot(x, agg["frac_all_three"], color="#1F4E79", linewidth=1.6, label="All three overlap")
        ax.plot(x, agg["frac_discordant"], color="#C43C39", linewidth=1.3, label="Discordant overlap")
        ax.plot(x, agg["frac_any"], color="#2E8B57", linewidth=1.0, linestyle="--", label="Any covered")
        ax.set_ylim(0, 1)
        ax.set_xlabel(f"{chrom} position (Mb)")
        ax.set_ylabel("Fraction of intervals")
        ax.set_title(f"Cohort overlap profile by bin ({chrom})")
        ax.grid(alpha=0.2, linewidth=0.5)
        ax.legend(loc="upper right", fontsize=8, frameon=False)
        plt.tight_layout()

        base = os.path.join(outdir, f"mosdepth_overlap_{safe_name(chrom)}_cohort_profile")
        out_png = f"{base}.png"
        out_svg = f"{base}.svg"
        fig.savefig(out_png, dpi=220)
        fig.savefig(out_svg)
        plt.close(fig)
        produced.extend([out_png, out_svg])
    return produced


def plot_chromosome_heatmaps(bin_df: pd.DataFrame, outdir: str, region_bin_size: int) -> List[str]:
    produced: List[str] = []
    chroms = sorted(bin_df["chrom"].dropna().unique(), key=chrom_sort_key)
    metric_defs = [
        ("frac_all_three", "all_three_overlap", "viridis_r", "Lower is worse overlap"),
        ("frac_discordant", "discordant_overlap", "magma", "Higher is worse overlap"),
    ]

    for chrom in chroms:
        sub = bin_df[bin_df["chrom"] == chrom].copy()
        if sub.empty:
            continue
        patients = sorted(sub["patient"].dropna().unique())
        bins_sorted = sorted(sub["bin_start"].dropna().unique())
        if not patients or not bins_sorted:
            continue

        for metric, metric_slug, cmap, subtitle in metric_defs:
            piv = (
                sub.pivot_table(index="patient", columns="bin_start", values=metric, aggfunc="mean")
                .reindex(index=patients, columns=bins_sorted)
                .copy()
            )
            if piv.empty:
                continue

            h = max(4.5, 0.28 * len(patients) + 1.8)
            w = max(10.0, min(20.0, 0.06 * len(bins_sorted) + 6.0))
            fig, ax = plt.subplots(figsize=(w, h))
            im = ax.imshow(piv.to_numpy(dtype=float), aspect="auto", vmin=0, vmax=1, cmap=cmap)
            ax.set_yticks(np.arange(len(patients)))
            ax.set_yticklabels(patients, fontsize=7)

            n_bins = len(bins_sorted)
            if n_bins <= 16:
                tick_idx = np.arange(n_bins)
            else:
                tick_idx = np.linspace(0, n_bins - 1, num=16, dtype=int)
            tick_bins = [bins_sorted[i] for i in tick_idx]
            tick_labels = [f"{(b + (region_bin_size / 2.0)) / 1e6:.1f}" for b in tick_bins]
            ax.set_xticks(tick_idx)
            ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=7)
            ax.set_xlabel(f"{chrom} position (Mb)")
            ax.set_title(f"{chrom}: {metric_slug} by patient and bin ({subtitle})", fontsize=10)

            cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
            cbar.set_label("Fraction (0-1)")
            plt.tight_layout()

            base = os.path.join(outdir, f"mosdepth_overlap_{safe_name(chrom)}_{metric_slug}_heatmap")
            out_png = f"{base}.png"
            out_svg = f"{base}.svg"
            fig.savefig(out_png, dpi=220)
            fig.savefig(out_svg)
            plt.close(fig)
            produced.extend([out_png, out_svg])
    return produced


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot overlap between DNA normal, DNA tumor, RNA tumor using mosdepth md.regions.bed.gz")
    ap.add_argument("--samples", required=True, help="SAMPLES.tsv path")
    ap.add_argument("--mosdepth-dir", required=True, help="Directory containing mosdepth patient folders")
    ap.add_argument("--dna-normal-label", default="DNA_NORMAL")
    ap.add_argument("--dna-tumor-label", default="DNA_TUMOR")
    ap.add_argument("--rna-tumor-label", default="RNA_TUMOR")
    ap.add_argument("--depth-threshold", type=float, default=10.0, help="Coverage threshold for overlap")
    ap.add_argument("--region-bin-size", type=int, default=250000, help="Genomic bin size for chromosome-level overlap plots")
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
    bin_rows: List[pd.DataFrame] = []
    missing_rows: List[Dict[str, str]] = []

    for patient in patients:
        try:
            row, by_cat, by_bin = compute_patient_overlap(
                patient=patient,
                mosdepth_dir=args.mosdepth_dir,
                dna_normal_label=args.dna_normal_label,
                dna_tumor_label=args.dna_tumor_label,
                rna_tumor_label=args.rna_tumor_label,
                depth_threshold=args.depth_threshold,
                region_bin_size=args.region_bin_size,
            )
            summary_rows.append(row)
            category_rows.append(by_cat)
            bin_rows.append(by_bin)
            print(f"[ok] {patient}")
        except FileNotFoundError as e:
            missing_rows.append({"patient": patient, "missing_file": str(e)})
            print(f"[skip] {patient}: missing {e}")
        except ValueError as e:
            missing_rows.append({"patient": patient, "missing_file": str(e)})
            print(f"[skip] {patient}: {e}")

    if not summary_rows:
        raise SystemExit("ERROR: no patients processed successfully (all missing mosdepth inputs?)")

    summary_df = pd.DataFrame(summary_rows).sort_values("patient").reset_index(drop=True)
    categories_df = pd.concat(category_rows, ignore_index=True)
    bins_df = pd.concat(bin_rows, ignore_index=True)

    summary_tsv = os.path.join(args.outdir, "mosdepth_overlap_summary.tsv")
    categories_tsv = os.path.join(args.outdir, "mosdepth_overlap_categories.tsv")
    bins_tsv = os.path.join(args.outdir, "mosdepth_overlap_bins.tsv")
    worst_bins_tsv = os.path.join(args.outdir, "mosdepth_overlap_worst_bins.tsv")
    missing_tsv = os.path.join(args.outdir, "mosdepth_overlap_missing_inputs.tsv")
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    categories_df.to_csv(categories_tsv, sep="\t", index=False)
    bins_df.to_csv(bins_tsv, sep="\t", index=False)
    worst_bins_df = bins_df.sort_values(
        ["frac_discordant", "frac_all_three", "patient", "chrom", "bin_start"],
        ascending=[False, True, True, True, True],
    ).copy()
    worst_bins_df.to_csv(worst_bins_tsv, sep="\t", index=False)
    pd.DataFrame(missing_rows).to_csv(missing_tsv, sep="\t", index=False)

    stacked_png = os.path.join(args.outdir, "mosdepth_overlap_stacked.png")
    stacked_svg = os.path.join(args.outdir, "mosdepth_overlap_stacked.svg")
    jaccard_png = os.path.join(args.outdir, "mosdepth_overlap_jaccard_heatmap.png")
    jaccard_svg = os.path.join(args.outdir, "mosdepth_overlap_jaccard_heatmap.svg")

    plot_stacked_overlap(categories_df, stacked_png, stacked_svg)
    plot_jaccard_heatmap(summary_df, jaccard_png, jaccard_svg)
    removed_plots = clear_chromosome_plot_outputs(args.outdir)
    chrom_profile_files = plot_chromosome_profiles(bins_df, args.outdir, args.region_bin_size)
    chrom_heatmap_files = plot_chromosome_heatmaps(bins_df, args.outdir, args.region_bin_size)

    print(f"[done] summary: {summary_tsv}")
    print(f"[done] categories: {categories_tsv}")
    print(f"[done] bins: {bins_tsv}")
    print(f"[done] worst bins: {worst_bins_tsv}")
    print(f"[done] missing: {missing_tsv}")
    print(f"[done] plots: {stacked_png}, {stacked_svg}, {jaccard_png}, {jaccard_svg}")
    print(f"[done] removed stale chromosome plot files: {removed_plots}")
    print(f"[done] chromosome cohort profiles: {len(chrom_profile_files)} files")
    print(f"[done] chromosome heatmaps: {len(chrom_heatmap_files)} files")


if __name__ == "__main__":
    main()
