#!/usr/bin/env python3
import argparse
import csv
import os
from collections import OrderedDict


def load_mupexi_rows(path, event_type):
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
    return header or [], rows


def write_tsv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def gather(inputs, event_type):
    rows = []
    fieldnames = []
    seen = OrderedDict()
    for patient, path in inputs:
        header, file_rows = load_mupexi_rows(path, event_type)
        if header and not fieldnames:
            fieldnames = ['patient_id', 'event_type', 'event_id', 'source_set', 'source_file'] + header
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
    return fieldnames, rows


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
    ap.add_argument('--snv-outfile', required=True)
    ap.add_argument('--fus-outfile', required=True)
    args = ap.parse_args()

    snv_inputs = parse_map(args.snv_input)
    fus_inputs = parse_map(args.fus_input)

    snv_fields, snv_rows = gather(snv_inputs, 'SNV')
    fus_fields, fus_rows = gather(fus_inputs, 'FUS')

    write_tsv(args.snv_outfile, snv_rows, snv_fields)
    write_tsv(args.fus_outfile, fus_rows, fus_fields)

    print(f'[done] gathered SNV rows={len(snv_rows)} -> {args.snv_outfile}')
    print(f'[done] gathered FUS rows={len(fus_rows)} -> {args.fus_outfile}')


if __name__ == '__main__':
    main()
