#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash research/extract_rna_editing_clusters_from_config.sh -c CONFIG [-s PATIENT] [--max-distance N] [--min-cluster-size N] [-o OUTDIR]
USAGE
}

outdir="research/output"
sample=""
max_distance=500
min_cluster_size=2

while :; do
  case ${1:-} in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --max-distance) max_distance="$2"; shift 2 ;;
    --min-cluster-size) min_cluster_size="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
done

[ -n "${config:-}" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
phased_ext="${rna7_phased_vcf_extension:-${phased_vcf_extension:-}}"
[ -n "$phased_ext" ] || { echo "ERROR: phased extension not found in CONFIG (rna7_phased_vcf_extension or phased_vcf_extension)" >&2; exit 1; }

mkdir -p "$outdir"
inputs=()

declare -A seen=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac
  sid=$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')
  patient="$sid"
  for lab in "${dna_normal_label:-DNA_NORMAL}" "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}" "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR"; do
    patient="${patient%_${lab}}"
  done
  [ -n "$patient" ] || continue
  [ -n "${seen[$patient]:-}" ] && continue
  seen["$patient"]=1

  if [ -n "$sample" ] && [ "$sample" != "$sid" ] && [ "$sample" != "$patient" ]; then
    continue
  fi

  vcf="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}/${patient}_${phased_ext}"
  if [ -f "$vcf" ]; then
    inputs+=("--input" "${patient}=${vcf}")
  else
    echo "[warn] missing phased VCF, skip: $vcf" >&2
  fi
done < "$samples"

if [ "${#inputs[@]}" -eq 0 ]; then
  echo "ERROR: no phased VCFs found" >&2
  exit 1
fi

python3 "$(dirname "$0")/extract_rna_editing_clusters.py" \
  "${inputs[@]}" \
  --out-variants "${outdir}/rna_edit_variants.tsv" \
  --out-clusters "${outdir}/rna_edit_clusters.tsv" \
  --max-distance "$max_distance" \
  --min-cluster-size "$min_cluster_size" \
  --rna-label "${rna_tumor_label:-RNA_TUMOR}" \
  --tumor-label "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" \
  --normal-label "${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"

echo "[done] wrote: ${outdir}/rna_edit_variants.tsv and ${outdir}/rna_edit_clusters.tsv"
