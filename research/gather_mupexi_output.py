#!/usr/bin/env python3
import argparse
import csv
import os
import re
from collections import OrderedDict

from annotate_dedup_vep import load_vep, parse_uploaded_variation


def normalize_chrom(value):
    text = str(value).strip()
    if text.lower().startswith('chr'):
        return text[3:]
    return text


def normalize_position_text(value):
    text = str(value).strip()
    if text in ('', 'NA'):
        return ''
    text = text.strip("[]'\"")
    return text


def normalize_protein_position(value):
    text = normalize_position_text(value)
    if not text:
        return ''
    text = re.sub(r'^\[|\]$', '', text)
    text = text.replace("'", "").replace('"', '')
    parts = [p.strip() for p in text.split(',') if p.strip()]
    if len(parts) == 1:
        return parts[0]
    return text


def parse_location_position(location):
    text = str(location).strip()
    if ':' not in text:
        return ''
    return normalize_position_text(text.split(':', 1)[1])


def parse_uploaded_parts(uploaded):
    parsed = parse_uploaded_variation(uploaded)
    if parsed is None:
        return '', '', '', ''
    chrom, pos, ref, alt = parsed
    return chrom, str(pos), ref, alt


def build_vep_lookup(path):
    _meta, _header, rows = load_vep(path)
    exact = {}
    fallback = {}
    for row in rows:
        uploaded = str(row.get('Uploaded_variation', '')).strip()
        if not uploaded:
            continue
        feature = str(row.get('Feature', '')).strip()
        aa = str(row.get('Amino_acids', '')).strip()
        chrom = normalize_chrom(str(row.get('Location', '')).split(':', 1)[0])
        pos_text = parse_location_position(row.get('Location', ''))
        prot = normalize_protein_position(row.get('Protein_position', ''))
        rec = {
            'uploaded_variation': uploaded,
            'Reference_Allele': 'NA',
            'Tumor_Seq_Allele2': 'NA',
        }
        _chrom2, _pos2, ref, alt = parse_uploaded_parts(uploaded)
        if ref:
            rec['Reference_Allele'] = ref
        if alt:
            rec['Tumor_Seq_Allele2'] = alt
        exact[(feature, chrom, pos_text, prot, aa)] = rec
        fallback[(feature, chrom, pos_text, aa)] = rec
    return exact, fallback


def annotate_snv_rows_with_vep(rows, vep_path):
    if not vep_path or not os.path.isfile(vep_path):
        for row in rows:
            row['uploaded_variation'] = 'NA'
            row['Reference_Allele'] = 'NA'
            row['Tumor_Seq_Allele2'] = 'NA'
        return 0
    exact, fallback = build_vep_lookup(vep_path)
    matched = 0
    for row in rows:
        feature = str(row.get('Transcript_ID', '')).strip()
        chrom = normalize_chrom(row.get('Chr', ''))
        pos_text = normalize_position_text(row.get('Genomic_Position', ''))
        prot = normalize_protein_position(row.get('Protein_position', ''))
        aa = str(row.get('Amino_Acid_Change', '')).strip()
        rec = exact.get((feature, chrom, pos_text, prot, aa))
        if rec is None:
            rec = fallback.get((feature, chrom, pos_text, aa))
        if rec is None:
            row['uploaded_variation'] = 'NA'
            row['Reference_Allele'] = 'NA'
            row['Tumor_Seq_Allele2'] = 'NA'
            continue
        row.update(rec)
        matched += 1
    return matched


def load_mupexi_rows(path, event_type, vep_path=''):
    header = None
    rows = []
    with open(path, 'r') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if header is None:
                header = cols
                continue
            row = {header[i]: cols[i] if i < len(cols) else '' for i in range(len(header))}
            if event_type == 'SNV':
                event_id = row.get('mutation_id_vep') or row.get('unique_peptide_id', '').split('|')[0]
            else:
                event_id = row.get('fusion_id', '')
            row['event_id'] = event_id or 'NA'
            rows.append(row)
    matched = 0
    if event_type == 'SNV':
        matched = annotate_snv_rows_with_vep(rows, vep_path)
        for row in rows:
            row['event_id'] = row.get('uploaded_variation') or row.get('event_id') or 'NA'
    return header or [], rows, matched


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def gather(inputs, event_type, vep_inputs=None):
    rows = []
    fieldnames = []
    seen = OrderedDict()
    total_matched = 0
    for patient, path in inputs:
        vep_path = ''
        if vep_inputs:
            vep_path = vep_inputs.get(patient, '')
        header, file_rows, matched = load_mupexi_rows(path, event_type, vep_path)
        total_matched += matched
        if header and not fieldnames:
            prefix_fields = ['patient_id', 'event_type', 'event_id', 'source_set', 'source_file']
            if event_type == 'SNV':
                prefix_fields += ['uploaded_variation', 'Reference_Allele', 'Tumor_Seq_Allele2']
            fieldnames = prefix_fields + header
        for row in file_rows:
            row = dict(row)
            row['patient_id'] = patient
            row['event_type'] = event_type
            row['source_set'] = row.get('Mutation_Origin', 'FUSION' if event_type == 'FUS' else 'NA')
            row['source_file'] = path
            key = (patient, row.get('event_id', ''), row.get('HLA_allele', ''), row.get('Mut_peptide', ''), row.get('Norm_peptide', ''))
            seen[key] = row
    rows = list(seen.values())
    rows.sort(key=lambda r: (r.get('patient_id', ''), r.get('event_id', ''), r.get('HLA_allele', ''), r.get('Mut_peptide', '')))
    if not fieldnames:
        fieldnames = ['patient_id', 'event_type', 'event_id', 'source_set', 'source_file']
    return fieldnames, rows, total_matched


def parse_map(items):
    out = []
    for item in items:
        if '=' not in item:
            raise SystemExit(f'ERROR: expected PATIENT=/path form, got: {item}')
        patient, path = item.split('=', 1)
        out.append((patient, path))
    return out


def main():
    ap = argparse.ArgumentParser(description='Gather MuPeXI SNV and fusion outputs into cohort TSVs.')
    ap.add_argument('--snv-input', action='append', default=[], help='PATIENT=/path/to/patient_snv.mupexi')
    ap.add_argument('--fus-input', action='append', default=[], help='PATIENT=/path/to/patient_fus.mupexi')
    ap.add_argument('--vep-input', action='append', default=[], help='PATIENT=/path/to/patient_vep.vep(.gz)')
    ap.add_argument('--snv-outfile', required=True)
    ap.add_argument('--fus-outfile', required=True)
    args = ap.parse_args()

    snv_inputs = parse_map(args.snv_input)
    fus_inputs = parse_map(args.fus_input)
    vep_inputs = dict(parse_map(args.vep_input))

    snv_fields, snv_rows, snv_matched = gather(snv_inputs, 'SNV', vep_inputs=vep_inputs)
    fus_fields, fus_rows, _fus_matched = gather(fus_inputs, 'FUS')

    write_tsv(args.snv_outfile, snv_rows, snv_fields)
    write_tsv(args.fus_outfile, fus_rows, fus_fields)

    print(f'[done] gathered SNV rows={len(snv_rows)} -> {args.snv_outfile}')
    print(f'[info] SNV rows with recovered genomic event_id from VEP={snv_matched}')
    print(f'[done] gathered FUS rows={len(fus_rows)} -> {args.fus_outfile}')


if __name__ == '__main__':
    main()
