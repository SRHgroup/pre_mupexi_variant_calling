#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash bin/check_outputs.sh -c CONFIG [-s SAMPLE]
USAGE
}

while :; do
  case ${1:-} in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c requires a path" >&2; exit 1; }
      config=$2
      shift
      ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s requires a sample" >&2; exit 1; }
      selected_sample=$2
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *) break ;;
  esac
  shift
done

[ -n "${config:-}" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }
source "$config"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label:-DNA_NORMAL}" "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}"
    "${out_dna_normal_label:-DNA_NORMAL}" "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do value="${value%_${label}}"; done
  printf '%s\n' "$value"
}

sample_is_requested() {
  local sample_id="$1"
  local patient_id="$2"
  local requested="${selected_sample:-}"
  [ -z "$requested" ] || [ "$sample_id" = "$requested" ] || [ "$patient_id" = "$requested" ]
}

missing=0
printf "sample\tstatus\tmissing_path\n"

declare -A seen_patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1

  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  outdir="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  expected=(
    "${vcfdir}/${name}_${output_extension_20}"
    "${vcfdir}/${name}_${output_extension_201}"
    "${vcfdir}/${name}_${output_extension_202}"
    "${vcfdir}/${name}_${output_extension_30}"
    "${outdir}/${name}_${rna_only_vcf_extension}"
    "${outdir}/${name}_${filtered_edit_labeled_vcf_extension}"
    "${outdir}/${name}_${annot_vcf_extension}"
    "${outdir}/${name}_${rna_summarised_vcf_extension}"
    "${outdir}/${name}_${rna_vcf_knownsites_extension}"
    "${outdir}/${name}_${merged_vcf_extension}"
    "${outdir}/${name}_${phased_vcf_extension}"
  )

  for out in "${expected[@]}"; do
    if [ ! -f "$out" ]; then
      printf "%s\tMISSING\t%s\n" "$name" "$out"
      missing=1
    fi
  done
done < "$samples"

if [ "$missing" -eq 0 ]; then
  echo "All expected files found."
fi
