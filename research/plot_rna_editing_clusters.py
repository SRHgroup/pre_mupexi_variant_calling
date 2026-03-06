#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def is_canonical(chrom):
    s = str(chrom).replace('chr', '')
    if s in ('X', 'Y'):
        return True
    return s.isdigit() and 1 <= int(s) <= 22


def expanded_copy_rows(df):
    rows = []
    for _, r in df.iterrows():
        gt = str(r.get('rna_gt', ''))
        ps = str(r.get('ps', ''))
        copies = []
        phased = False
        if '|' in gt:
            parts = gt.split('|')
            if len(parts) >= 2:
                # Alt on copy i if allele is non-reference and not missing.
                if parts[0] not in ('0', '.', ''):
                    copies.append(1)
                if parts[1] not in ('0', '.', ''):
                    copies.append(2)
                if not copies:
                    copies = [1, 2]
                # Phase state comes from GT separator. PS is only for block grouping.
                phased = True
            else:
                copies = ['U']
        else:
            # Unphased or missing genotype: show once in middle lane.
            copies = ['U']
        for cp in copies:
            rr = r.copy()
            rr['copy'] = cp
            rr['is_phased'] = phased
            rows.append(rr)
    return pd.DataFrame(rows)


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


def write_empty_plot(out_prefix, reason):
    fig, ax = plt.subplots(1, 1, figsize=(12, 3))
    ax.axis('off')
    ax.text(0.5, 0.5, f'No variants to plot\n{reason}', ha='center', va='center', fontsize=12)
    out_png = f"{out_prefix}.png"
    out_pdf = f"{out_prefix}.pdf"
    os.makedirs(os.path.dirname(out_png) or '.', exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.savefig(out_pdf)
    print(f"[warn] empty plot written: {out_png}, {out_pdf} ({reason}; check extract filters like --min-alt-count)")


def main():
    ap = argparse.ArgumentParser(description='Plot multi-patient RNA-editing clusters.')
    ap.add_argument('--variants', required=True)
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--top-clusters', type=int, default=20)
    ap.add_argument('--canonical-only', action='store_true', default=True)
    args = ap.parse_args()

    v = pd.read_csv(args.variants, sep='\t')
    c = pd.read_csv(args.clusters, sep='\t')
    if v.empty:
        write_empty_plot(args.out_prefix, 'variants TSV has no rows')
        return

    if args.canonical_only:
        v = v[v['chrom'].map(is_canonical)].copy()
        c = c[c['chrom'].map(is_canonical)].copy()
        if v.empty:
            write_empty_plot(args.out_prefix, 'no canonical-chromosome variants')
            return

    v['pos'] = pd.to_numeric(v['pos'])
    v['rna_dp'] = pd.to_numeric(v['rna_dp'], errors='coerce')
    chr_order = build_chr_order(v['chrom'].dropna().unique())
    patients = sorted(v['patient'].dropna().unique())
    p_rank = {p: i for i, p in enumerate(patients)}

    # Numeric pseudo-genomic axis by chr blocks.
    offsets = {}
    bounds = []
    cur = 0
    chr_gap = 20_000_000
    for ch in chr_order:
        max_pos = v.loc[v['chrom'] == ch, 'pos'].max()
        start = cur
        offsets[ch] = cur
        cur += int(max_pos) + chr_gap
        bounds.append((ch, start, cur))

    v['gpos'] = v.apply(lambda r: offsets.get(r['chrom'], 0) + r['pos'], axis=1)
    v['known_db_hit'] = v['known_db'].fillna('').astype(str).map(lambda x: x not in ('', '.'))

    # Map each variant to a cluster_id by interval overlap (patient/chrom/start-end).
    v['cluster_id'] = ''
    if not c.empty:
        c_idx = c[['cluster_id', 'patient', 'chrom', 'start', 'end']].copy()
        c_idx['start'] = pd.to_numeric(c_idx['start'], errors='coerce')
        c_idx['end'] = pd.to_numeric(c_idx['end'], errors='coerce')
        c_idx = c_idx.dropna(subset=['start', 'end'])
        for (pt, ch), g in c_idx.groupby(['patient', 'chrom']):
            mask = (v['patient'] == pt) & (v['chrom'] == ch)
            if not mask.any():
                continue
            iv = g.sort_values(['start', 'end'])
            pos = v.loc[mask, 'pos']
            cid = pd.Series('', index=pos.index, dtype=object)
            for _, r in iv.iterrows():
                hit = (pos >= int(r['start'])) & (pos <= int(r['end'])) & (cid == '')
                cid.loc[hit] = str(r['cluster_id'])
            v.loc[mask, 'cluster_id'] = cid

    # Build cluster color map with high contrast: adjacent clusters get distant palette colors.
    cluster_color_map = {}
    if not c.empty:
        c_work = c.copy()
        c_work['gstart'] = c_work.apply(lambda r: offsets.get(r['chrom'], 0) + int(r['start']), axis=1)
        c_work = c_work.sort_values(['patient', 'gstart', 'end'])
        palette = (
            list(plt.get_cmap('tab20').colors) +
            list(plt.get_cmap('Dark2').colors) +
            list(plt.get_cmap('Set1').colors) +
            list(plt.get_cmap('Accent').colors)
        )
        # Deduplicate identical tuples while preserving order.
        seen_colors = set()
        palette_unique = []
        for col in palette:
            if col not in seen_colors:
                palette_unique.append(col)
                seen_colors.add(col)
        palette = palette_unique
        # Jump through palette with a coprime step to maximize contrast between nearby clusters.
        step = 7
        ncol = len(palette)
        for i, cid in enumerate(c_work['cluster_id'].astype(str)):
            cluster_color_map[cid] = palette[(i * step) % ncol]

    vx = expanded_copy_rows(v)
    vx['patient_idx'] = vx['patient'].map(p_rank)
    copy_track_gap = 3
    vx['copy_row'] = vx['patient_idx'] * copy_track_gap + vx['copy'].map({1: 0.0, 'U': 0.5, 2: 1.0})
    copy_x_offset = 600000
    vx['gpos_copy'] = vx['gpos'] + vx['copy'].map({1: -copy_x_offset, 'U': 0, 2: copy_x_offset})
    # Small deterministic y-jitter to reduce overplotting while staying on each lane.
    rng = np.random.default_rng(42)
    vx['copy_row_plot'] = vx['copy_row'] + rng.uniform(-0.08, 0.08, size=len(vx))

    signature_marker = {'ADAR': 's', 'APOBEC3': 'o', 'OTHER': '^'}

    fig_h = max(5.5, min(10.5, 4.0 + 0.14 * len(patients)))
    fig, (ax_top, ax) = plt.subplots(2, 1, figsize=(16, fig_h), gridspec_kw={'height_ratios': [1, 3]}, sharex=True)

    # Top: depth density-ish trend via binned mean DP.
    bins = 300
    v_nonan = v.dropna(subset=['rna_dp']).copy()
    if not v_nonan.empty:
        v_nonan['bin'] = pd.cut(v_nonan['gpos'], bins=bins, labels=False)
        b = v_nonan.groupby('bin', as_index=False).agg(gpos=('gpos', 'mean'), mean_dp=('rna_dp', 'mean'))
        ax_top.plot(b['gpos'], b['mean_dp'], color='#2ca02c', linewidth=1.5)
        ax_top.fill_between(b['gpos'], b['mean_dp'], alpha=0.2, color='#2ca02c')
    ax_top.set_ylabel('Mean RNA DP')
    ax_top.set_title('RNA-editing Variant Landscape Across Patients', pad=30)

    # Bottom: signature by shape, cluster by color, known-db by alpha.
    for sig, sub in vx.groupby('signature'):
        marker = signature_marker.get(sig, 'o')
        phased = sub[sub['is_phased']]
        unphased = sub[~sub['is_phased']]
        for part, hollow in ((phased, False), (unphased, True)):
            if part.empty:
                continue
            for cid, g in part.groupby(part['cluster_id'].fillna('')):
                sub_known = g[g['known_db_hit']]
                sub_novel = g[~g['known_db_hit']]
                cval = cluster_color_map.get(str(cid), (0.72, 0.72, 0.72))
                if not sub_novel.empty:
                    ax.scatter(
                        sub_novel['gpos_copy'], sub_novel['copy_row_plot'],
                        s=9,
                        alpha=0.15,
                        marker=marker,
                        facecolors='none' if hollow else cval,
                        edgecolors=cval if hollow else 'black',
                        linewidths=0.7 if hollow else 0.25
                    )
                if not sub_known.empty:
                    ax.scatter(
                        sub_known['gpos_copy'], sub_known['copy_row_plot'],
                        s=9,
                        alpha=1.0,
                        marker=marker,
                        facecolors='none' if hollow else cval,
                        edgecolors=cval if hollow else 'black',
                        linewidths=0.7 if hollow else 0.25
                    )

    # Draw PS connectors to make haplotype blocks visually explicit.
    psub = vx[(vx['is_phased']) & (vx['ps'].notna()) & (vx['ps'].astype(str) != '') & (vx['ps'].astype(str) != '.')]
    if not psub.empty:
        for (_, _, _, _), g in psub.groupby(['patient', 'chrom', 'copy', 'ps']):
            if len(g) < 2:
                continue
            g = g.sort_values('gpos_copy')
            ax.plot(g['gpos_copy'], g['copy_row'], color='#444444', linewidth=0.8, alpha=0.35)

    # Highlight top clusters by size.
    if not c.empty:
        c2 = c.sort_values('n_variants', ascending=False).head(args.top_clusters)
        for _, r in c2.iterrows():
            x0 = offsets.get(r['chrom'], 0) + int(r['start'])
            x1 = offsets.get(r['chrom'], 0) + int(r['end'])
            pidx = p_rank.get(r['patient'], None)
            if pidx is None:
                continue
            y1 = pidx * copy_track_gap
            y2 = pidx * copy_track_gap + 1
            ax.plot([x0 - copy_x_offset, x1 - copy_x_offset], [y1, y1], color='black', linewidth=1.8, alpha=0.35)
            ax.plot([x0 + copy_x_offset, x1 + copy_x_offset], [y2, y2], color='black', linewidth=1.8, alpha=0.35)

    # Visual chromosome separation.
    for i, (_, start, end) in enumerate(bounds):
        if i % 2 == 1:
            ax.axvspan(start, end, color='#f3f3f3', alpha=0.35, zorder=0)
            ax_top.axvspan(start, end, color='#f3f3f3', alpha=0.35, zorder=0)
        if i > 0:
            ax.axvline(start, color='black', alpha=0.18, linewidth=0.8)
            ax_top.axvline(start, color='black', alpha=0.18, linewidth=0.8)

    yticks = []
    ylabels = []
    for i, p in enumerate(patients):
        yticks.extend([i * copy_track_gap, i * copy_track_gap + 1])
        ylabels.extend([f'{p} copy1', f'{p} copy2'])
        # subtle separator between patients
        if i > 0:
            ax.axhline(i * copy_track_gap - 1, color='black', alpha=0.08, linewidth=0.8)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_ylim(-0.8, max(1.8, len(patients) * copy_track_gap - 0.2))
    ax.set_ylabel('Patient / Copy')
    ax.set_xlabel('Genomic Position (chromosome-concatenated)')
    sig_handles = [
        Line2D([], [], marker='s', linestyle='None', color='black', label='ADAR', markersize=5),
        Line2D([], [], marker='o', linestyle='None', color='black', label='APOBEC3', markersize=5),
        Line2D([], [], marker='^', linestyle='None', color='black', label='OTHER', markersize=5),
    ]
    shape_handles = [
        Line2D([], [], marker='o', linestyle='None', color='black', alpha=1.0, label='KNOWN_RNAEDIT_DB hit', markersize=5),
        Line2D([], [], marker='o', linestyle='None', color='black', alpha=0.15, label='No KNOWN_RNAEDIT_DB', markersize=5),
    ]
    phase_handles = [
        Line2D([], [], marker='o', linestyle='None', color='black', label='Phased', markersize=4),
        Line2D([], [], marker='o', linestyle='None', markerfacecolor='none', color='black', label='Unphased (middle lane)', markersize=4),
    ]
    legend_y = 1.22
    leg1 = ax_top.legend(
        handles=sig_handles,
        loc='upper left',
        bbox_to_anchor=(0.0, legend_y),
        ncol=3,
        frameon=False,
        title='Signature'
    )
    ax_top.add_artist(leg1)
    leg2 = ax_top.legend(
        handles=shape_handles,
        loc='upper center',
        bbox_to_anchor=(0.5, legend_y),
        frameon=False,
        title='Database'
    )
    ax_top.add_artist(leg2)
    ax_top.legend(
        handles=phase_handles,
        loc='upper right',
        bbox_to_anchor=(1.0, legend_y),
        frameon=False,
        title='Phasing'
    )

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

    plt.tight_layout(rect=(0, 0, 1, 0.92))
    out_png = f"{args.out_prefix}.png"
    out_pdf = f"{args.out_prefix}.pdf"
    os.makedirs(os.path.dirname(out_png) or '.', exist_ok=True)
    plt.savefig(out_png, dpi=180)
    plt.savefig(out_pdf)
    print(f"[done] plots: {out_png}, {out_pdf}")


if __name__ == '__main__':
    main()
