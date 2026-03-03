#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  rna_minus_dna_by_pos.sh -r RNA.vcf.gz -d DNA.vcf.gz [-o OUT.vcf.gz] [-t THREADS]

What it does (ONLY):
  Removes every variant from RNA whose CHROM+POS exists in DNA (ignores REF/ALT).
  No PASS filtering. No normalization. No QC. No indexing.

Required:
  -r, --rna     RNA VCF (.vcf.gz)
  -d, --dna     DNA VCF (.vcf.gz)

Optional:
  -o, --out     Output VCF (.vcf.gz). Default: <RNA>.minusDNApos.vcf.gz next to RNA
  -t, --threads Threads for bcftools (default: 4)
  -h, --help    Show help

Notes:
  - If OUT directory doesn't exist, the script will fail (it will not create dirs).
  - Contig naming must match between files (chr1 vs 1).
EOF
}

die(){ echo "ERROR: $*" >&2; exit 1; }

RNA=""
DNA=""
OUT=""
THREADS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--rna) RNA="$2"; shift 2;;
    -d|--dna) DNA="$2"; shift 2;;
    -o|--out) OUT="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1 (use -h)";;
  esac
done

[[ -n "$RNA" ]] || die "Missing --rna"
[[ -n "$DNA" ]] || die "Missing --dna"
[[ -f "$RNA" ]] || die "RNA not found: $RNA"
[[ -f "$DNA" ]] || die "DNA not found: $DNA"
command -v bcftools >/dev/null 2>&1 || die "bcftools not in PATH"
[[ "$RNA" == *.vcf.gz ]] || die "RNA must be .vcf.gz"
[[ "$DNA" == *.vcf.gz ]] || die "DNA must be .vcf.gz"

if [[ -z "$OUT" ]]; then
  OUT="${RNA%.vcf.gz}.minusDNApos.vcf.gz"
fi
[[ "$OUT" == *.vcf.gz ]] || die "--out must end with .vcf.gz"

TMPBED="$(mktemp "${TMPDIR:-/tmp}/dna_pos.XXXXXX.bed")"
trap 'rm -f "$TMPBED"' EXIT

# Build 0-based BED intervals covering exactly each POS
bcftools query -f '%CHROM\t%POS\n' "$DNA" \
| awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2}' \
| sort -u > "$TMPBED"

# Exclude those positions from RNA
bcftools view --threads "$THREADS" -T "^$TMPBED" -Oz -o "$OUT" "$RNA"

echo "OK: wrote $OUT"
