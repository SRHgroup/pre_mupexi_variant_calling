#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
from collections import defaultdict, Counter


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_info(info):
    out = {}
    if info == '.' or info == '':
        return out
    for token in info.split(';'):
        if '=' in token:
            k, v = token.split('=', 1)
            out[k] = v
        else:
            out[token] = True
    return out


def detect_label_sample(samples, requested, labels):
    if requested:
        for s in samples:
            if s == requested:
                return s
    patterns = []
    for lab in labels:
        patterns.append(re.compile(rf"(?:_|^){re.escape(lab)}(?:\.|$)", re.IGNORECASE))
        if 'TUMOR' in lab:
            alt = lab.replace('TUMOR', 'TUMOUR')
            patterns.append(re.compile(rf"(?:_|^){re.escape(alt)}(?:\.|$)", re.IGNORECASE))
        if 'TUMOUR' in lab:
            alt = lab.replace('TUMOUR', 'TUMOR')
            patterns.append(re.compile(rf"(?:_|^){re.escape(alt)}(?:\.|$)", re.IGNORECASE))
    for s in samples:
        for p in patterns:
            if p.search(s):
                return s
    return None


def detect_non_normal_sample(samples, normal_labels):
    pats = []
    for lab in normal_labels:
        pats.append(re.compile(rf"(?:_|^){re.escape(lab)}(?:\\.|$)", re.IGNORECASE))
    for s in samples:
        is_norm = any(p.search(s) for p in pats)
        if not is_norm:
            return s
    return None


def parse_format_value(fmt_keys, sample_value, key):
    vals = sample_value.split(':')
    idx = {k: i for i, k in enumerate(fmt_keys)}
    if key not in idx:
        return None
    i = idx[key]
    if i >= len(vals):
        return None
    v = vals[i]
    if v in ('.', ''):
        return None
    return v


def parse_numeric(x):
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        return None


def infer_af_from_ad(ad):
    if not ad:
        return None
    parts = ad.split(',')
    if len(parts) < 2:
        return None
    try:
        ref = float(parts[0])
        alt = float(parts[1])
    except Exception:
        return None
    denom = ref + alt
    if denom <= 0:
        return None
    return alt / denom


def infer_alt_count_from_fields(ad, af, dp):
    if ad:
        parts = ad.split(',')
        if len(parts) >= 2:
            try:
                return float(parts[1])
            except Exception:
                pass
    if af is not None and dp is not None:
        try:
            return float(af) * float(dp)
        except Exception:
            return None
    return None


def infer_signature(info_map):
    if 'ADAR_SIG' in info_map:
        return 'ADAR'
    if 'APOBEC3_SIG' in info_map:
        return 'APOBEC3'
    return 'OTHER'


def infer_gene_and_transcript(info_map):
    # Try ANN= (SnpEff) first, then CSQ=.
    gene = ''
    transcript = ''
    if 'ANN' in info_map:
        ann = info_map['ANN'].split(',')[0]
        fields = ann.split('|')
        if len(fields) > 3:
            gene = fields[3]
        if len(fields) > 6:
            transcript = fields[6]
    elif 'CSQ' in info_map:
        csq = info_map['CSQ'].split(',')[0]
        fields = csq.split('|')
        # VEP layout varies; keep best-effort.
        for f in fields:
            if not transcript and re.match(r'^ENST\d+', f):
                transcript = f
            if not gene and re.match(r'^[A-Za-z0-9_.-]{2,}$', f) and not f.startswith('ENST'):
                gene = f
                break
    return gene, transcript


def passes_rna_edit_filter(info_map):
    if 'SOURCE_SET' in info_map and 'RNA_EDIT' in str(info_map['SOURCE_SET']):
        return True
    if 'ADAR_SIG' in info_map or 'APOBEC3_SIG' in info_map:
        return True
    if 'KNOWN_RNAEDIT_DB' in info_map and str(info_map['KNOWN_RNAEDIT_DB']) not in ('', '.'):
        return True
    return False


def build_clusters(variants, max_distance, min_size):
    clusters = []
    by_key = defaultdict(list)
    for v in variants:
        by_key[(v['patient'], v['chrom'])].append(v)

    cid = 0
    for (patient, chrom), rows in by_key.items():
        rows.sort(key=lambda r: r['pos'])
        cur = [rows[0]] if rows else []
        for r in rows[1:]:
            if r['pos'] - cur[-1]['pos'] <= max_distance:
                cur.append(r)
            else:
                if len(cur) >= min_size:
                    cid += 1
                    clusters.append(make_cluster_row(cid, patient, chrom, cur))
                cur = [r]
        if cur and len(cur) >= min_size:
            cid += 1
            clusters.append(make_cluster_row(cid, patient, chrom, cur))

    return clusters


def make_cluster_row(cid, patient, chrom, rows):
    positions = [r['pos'] for r in rows]
    sig = Counter(r['signature'] for r in rows)
    genes = Counter(r['gene'] for r in rows if r['gene'])
    tx = Counter(r['transcript'] for r in rows if r['transcript'])
    ps_vals = [r['ps'] for r in rows if r['ps']]
    phased = 1 if ps_vals and len(set(ps_vals)) == 1 and len(rows) >= 2 else 0
    mean_dp = sum(r['rna_dp'] for r in rows if r['rna_dp'] is not None) / max(1, sum(1 for r in rows if r['rna_dp'] is not None))

    return {
        'cluster_id': f'CL{cid:06d}',
        'patient': patient,
        'chrom': chrom,
        'start': min(positions),
        'end': max(positions),
        'n_variants': len(rows),
        'adar_count': sig.get('ADAR', 0),
        'apobec_count': sig.get('APOBEC3', 0),
        'other_count': sig.get('OTHER', 0),
        'mean_rna_dp': round(mean_dp, 3) if mean_dp == mean_dp else '',
        'single_ps_block': phased,
        'top_gene': genes.most_common(1)[0][0] if genes else '',
        'top_transcript': tx.most_common(1)[0][0] if tx else '',
    }


def main():
    ap = argparse.ArgumentParser(description='Extract RNA-editing variant clusters from phased VCFs.')
    ap.add_argument('--input', action='append', required=True, help='PATIENT=/path/to/phased.vcf.gz (repeatable)')
    ap.add_argument('--out-variants', required=True)
    ap.add_argument('--out-clusters', required=True)
    ap.add_argument('--max-distance', type=int, default=500)
    ap.add_argument('--min-cluster-size', type=int, default=2)
    ap.add_argument('--min-alt-count', type=float, default=10.0)
    ap.add_argument('--rna-sample', default='')
    ap.add_argument('--rna-label', action='append', default=['RNA_TUMOR'])
    ap.add_argument('--tumor-label', action='append', default=['TUMOR', 'DNA_TUMOR', 'DNA_TUMOUR'])
    ap.add_argument('--normal-label', action='append', default=['DNA_NORMAL'])
    args = ap.parse_args()

    pairs = []
    for item in args.input:
        if '=' not in item:
            raise SystemExit(f'--input must be PATIENT=VCF, got: {item}')
        p, v = item.split('=', 1)
        pairs.append((p, v))

    variants = []
    for patient, vcf in pairs:
        if not os.path.exists(vcf):
            print(f'[warn] missing input, skip: {patient} {vcf}')
            continue

        with open_text(vcf) as fh:
            header_samples = None
            sig_col_idx = None
            sig_sample_name = None
            for line in fh:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    cols = line.rstrip('\n').split('\t')
                    header_samples = cols[9:]
                    # Prefer RNA-labelled column. If absent (e.g., only DNA_NORMAL+TUMOR), use tumor/sample-not-normal column.
                    sig_sample_name = detect_label_sample(header_samples, args.rna_sample, args.rna_label)
                    if sig_sample_name is None:
                        sig_sample_name = detect_label_sample(header_samples, '', args.tumor_label)
                    if sig_sample_name is None:
                        sig_sample_name = detect_non_normal_sample(header_samples, args.normal_label)
                    if sig_sample_name is None:
                        print(f'[warn] no signal sample detected for {patient} in {vcf}')
                        break
                    sig_col_idx = 9 + header_samples.index(sig_sample_name)
                    continue
                if header_samples is None or sig_col_idx is None:
                    continue

                f = line.rstrip('\n').split('\t')
                if len(f) < 10:
                    continue
                chrom = f[0]
                pos = int(f[1])
                ref = f[3]
                alt = f[4]
                info_map = parse_info(f[7])

                if not passes_rna_edit_filter(info_map):
                    continue

                fmt_keys = f[8].split(':')
                sig_sv = f[sig_col_idx]
                gt = parse_format_value(fmt_keys, sig_sv, 'GT') or ''
                ad = parse_format_value(fmt_keys, sig_sv, 'AD') or ''
                af = parse_numeric(parse_format_value(fmt_keys, sig_sv, 'AF'))
                if af is None:
                    af = infer_af_from_ad(ad)
                dp = parse_numeric(parse_format_value(fmt_keys, sig_sv, 'DP'))
                alt_count = infer_alt_count_from_fields(ad, af, dp)
                if alt_count is None or alt_count < args.min_alt_count:
                    continue
                ps = parse_format_value(fmt_keys, sig_sv, 'PS') or ''
                pid = parse_format_value(fmt_keys, sig_sv, 'PID') or ''
                gene, transcript = infer_gene_and_transcript(info_map)

                variants.append({
                    'patient': patient,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'signature': infer_signature(info_map),
                    'rna_sample': sig_sample_name,
                    'rna_gt': gt,
                    'rna_ad': ad,
                    'rna_alt_count': round(alt_count, 3),
                    'rna_af': '' if af is None else round(af, 6),
                    'rna_dp': None if dp is None else float(dp),
                    'ps': ps,
                    'pid': pid,
                    'gene': gene,
                    'transcript': transcript,
                    'source_set': info_map.get('SOURCE_SET', ''),
                    'known_db': info_map.get('KNOWN_RNAEDIT_DB', ''),
                })

    clusters = build_clusters(variants, args.max_distance, args.min_cluster_size)

    os.makedirs(os.path.dirname(args.out_variants) or '.', exist_ok=True)
    os.makedirs(os.path.dirname(args.out_clusters) or '.', exist_ok=True)

    v_fields = [
        'patient', 'chrom', 'pos', 'ref', 'alt', 'signature', 'rna_sample', 'rna_gt', 'rna_ad', 'rna_alt_count', 'rna_af', 'rna_dp',
        'ps', 'pid', 'gene', 'transcript', 'source_set', 'known_db'
    ]
    with open(args.out_variants, 'w', newline='') as outv:
        w = csv.DictWriter(outv, fieldnames=v_fields, delimiter='\t')
        w.writeheader()
        for r in sorted(variants, key=lambda x: (x['patient'], x['chrom'], x['pos'])):
            rr = dict(r)
            if rr['rna_dp'] is None:
                rr['rna_dp'] = ''
            w.writerow(rr)

    c_fields = [
        'cluster_id', 'patient', 'chrom', 'start', 'end', 'n_variants', 'adar_count', 'apobec_count', 'other_count',
        'mean_rna_dp', 'single_ps_block', 'top_gene', 'top_transcript'
    ]
    with open(args.out_clusters, 'w', newline='') as outc:
        w = csv.DictWriter(outc, fieldnames=c_fields, delimiter='\t')
        w.writeheader()
        for c in sorted(clusters, key=lambda x: (x['patient'], x['chrom'], x['start'])):
            w.writerow(c)

    print(f'[done] variants={len(variants)} -> {args.out_variants}')
    print(f'[done] clusters={len(clusters)} -> {args.out_clusters}')


if __name__ == '__main__':
    main()
