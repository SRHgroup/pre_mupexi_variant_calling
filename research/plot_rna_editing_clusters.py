#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt


def build_chr_order(chrs):
    order = []
    for c in chrs:
        s = c.replace('chr', '')
        if s.isdigit():
            order.append((0, int(s), c))
        elif s == 'X':
            order.append((1, 23, c))
        elif s == 'Y':
            order.append((1, 24, c))
        else:
            order.append((2, 999, c))
    return [x[2] for x in sorted(order)]


def main():
    ap = argparse.ArgumentParser(description='Plot multi-patient RNA-editing clusters.')
    ap.add_argument('--variants', required=True)
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--top-clusters', type=int, default=20)
    args = ap.parse_args()

    v = pd.read_csv(args.variants, sep='\t')
    c = pd.read_csv(args.clusters, sep='\t')
    if v.empty:
        raise SystemExit('No variants to plot.')

    v['pos'] = pd.to_numeric(v['pos'])
    v['rna_dp'] = pd.to_numeric(v['rna_dp'], errors='coerce')
    chr_order = build_chr_order(v['chrom'].dropna().unique())
    chr_rank = {ch: i for i, ch in enumerate(chr_order)}
    patients = sorted(v['patient'].dropna().unique())
    p_rank = {p: i for i, p in enumerate(patients)}

    # Numeric pseudo-genomic axis by chr blocks.
    offsets = {}
    cur = 0
    for ch in chr_order:
        max_pos = v.loc[v['chrom'] == ch, 'pos'].max()
        offsets[ch] = cur
        cur += int(max_pos) + 1_000_000

    v['gpos'] = v.apply(lambda r: offsets.get(r['chrom'], 0) + r['pos'], axis=1)
    v['patient_idx'] = v['patient'].map(p_rank)

    color_map = {'ADAR': '#1f77b4', 'APOBEC3': '#d62728', 'OTHER': '#7f7f7f'}

    fig, (ax_top, ax) = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw={'height_ratios': [1, 4]}, sharex=True)

    # Top: depth density-ish trend via binned mean DP.
    bins = 300
    v_nonan = v.dropna(subset=['rna_dp']).copy()
    if not v_nonan.empty:
        v_nonan['bin'] = pd.cut(v_nonan['gpos'], bins=bins, labels=False)
        b = v_nonan.groupby('bin', as_index=False).agg(gpos=('gpos', 'mean'), mean_dp=('rna_dp', 'mean'))
        ax_top.plot(b['gpos'], b['mean_dp'], color='#2ca02c', linewidth=1.5)
        ax_top.fill_between(b['gpos'], b['mean_dp'], alpha=0.2, color='#2ca02c')
    ax_top.set_ylabel('Mean RNA DP')
    ax_top.set_title('RNA-editing Variant Landscape Across Patients')

    # Bottom: points by patient and signature.
    for sig, sub in v.groupby('signature'):
        ax.scatter(sub['gpos'], sub['patient_idx'], s=16, alpha=0.8, c=color_map.get(sig, '#333333'), label=sig)

    # Highlight top clusters by size.
    if not c.empty:
        c2 = c.sort_values('n_variants', ascending=False).head(args.top_clusters)
        for _, r in c2.iterrows():
            x0 = offsets.get(r['chrom'], 0) + int(r['start'])
            x1 = offsets.get(r['chrom'], 0) + int(r['end'])
            y = p_rank.get(r['patient'], None)
            if y is None:
                continue
            ax.plot([x0, x1], [y, y], color='black', linewidth=2, alpha=0.45)

    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels(patients)
    ax.set_ylabel('Patient')
    ax.set_xlabel('Genomic Position (chromosome-concatenated)')
    ax.legend(loc='upper right', ncol=3, frameon=False)

    # Chromosome tick labels
    xticks = []
    xlabels = []
    for ch in chr_order:
        vals = v.loc[v['chrom'] == ch, 'gpos']
        if vals.empty:
            continue
        xticks.append((vals.min() + vals.max()) / 2)
        xlabels.append(ch)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha='right')

    plt.tight_layout()
    out_png = f"{args.out_prefix}.png"
    out_pdf = f"{args.out_prefix}.pdf"
    os.makedirs(os.path.dirname(out_png) or '.', exist_ok=True)
    plt.savefig(out_png, dpi=180)
    plt.savefig(out_pdf)
    print(f"[done] plots: {out_png}, {out_pdf}")


if __name__ == '__main__':
    main()
