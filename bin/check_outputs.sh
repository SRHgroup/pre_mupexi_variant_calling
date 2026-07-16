#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash bin/check_outputs.sh -c CONFIG [-s SAMPLE] [-m MODE]
MODE: all (default), rna, germline, splicing
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
    -m|--mode)
      [ -n "${2:-}" ] || { echo "ERROR: -m requires a mode: all|rna|germline|splicing" >&2; exit 1; }
      mode=$2
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
mode="${mode:-all}"
case "$mode" in
  all|rna|germline|splicing) ;;
  *) echo "ERROR: invalid mode '$mode' (use all|rna|germline|splicing)" >&2; exit 1 ;;
esac
source "$config"

: "${samples:?CONFIG must define samples}"
if [ "$mode" = "splicing" ]; then
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  # shellcheck disable=SC1091
  source "${script_dir}/splicing_check_helpers.sh"
  splicing_check_star_root >/dev/null || {
    echo "ERROR: CONFIG must define splicing_sjdir, splicing_stardir, stardir, or bamdir" >&2
    exit 1
  }
  splicing_check_output_root >/dev/null || {
    echo "ERROR: CONFIG must define splicing_outdir or datadir" >&2
    exit 1
  }
else
  : "${vcfdir:?CONFIG must define vcfdir}"
  : "${bamdir:?CONFIG must define bamdir}"
  : "${gdna1_vcf_extension:?CONFIG must define gdna1_vcf_extension}"
  : "${gdna2_vcf_extension:?CONFIG must define gdna2_vcf_extension}"
  : "${gdna3_vcf_extension:?CONFIG must define gdna3_vcf_extension}"
  : "${gdna4_vcf_extension:?CONFIG must define gdna4_vcf_extension}"
  : "${rna1_vcf_extension:?CONFIG must define rna1_vcf_extension}"
  : "${rna2_labeled_vcf_extension:?CONFIG must define rna2_labeled_vcf_extension}"
  : "${rna3_knownsites_vcf_extension:?CONFIG must define rna3_knownsites_vcf_extension}"
  : "${rna4_summarised_vcf_extension:?CONFIG must define rna4_summarised_vcf_extension}"
  : "${rna5_qced_vcf_extension:?CONFIG must define rna5_qced_vcf_extension}"
  : "${rna6_merged_vcf_extension:?CONFIG must define rna6_merged_vcf_extension}"
  : "${rna7_smfixed_bam_suffix:?CONFIG must define rna7_smfixed_bam_suffix}"
  : "${rna7_phased_vcf_extension:?CONFIG must define rna7_phased_vcf_extension}"
fi

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

failed=0
checked_samples=0
splicing_steps=(spl1 spl2 spl3)
if [ -n "${splicing_normal_junction_ref:-}" ] || [[ "${splicing_spl4_input_suffix:-}" == *.spl3.5.* ]]; then
  splicing_steps+=(spl3.5)
fi
splicing_steps+=(spl4)
if [ "$mode" = "splicing" ]; then
  printf "sample\tstep\tstatus\tpath\n"
else
  printf "sample\tstatus\tpath\n"
fi

declare -A seen_patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  if [ "$mode" = "splicing" ]; then
    splicing_check_is_rna_sample_id "$sample_name" || continue
  fi
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1
  checked_samples=$((checked_samples + 1))

  if [ "$mode" = "splicing" ]; then
    for step in "${splicing_steps[@]}"; do
      completion="$(splicing_check_step_completion "$name" "$step")" || true
      status="${completion%%$'\t'*}"
      path="${completion#*$'\t'}"
      printf "%s\t%s\t%s\t%s\n" "$name" "$step" "$status" "$path"
      [ "$status" = "DONE" ] || failed=1
    done
    continue
  fi

  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  outdir="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"
  germline_dir="${vcfdir}/${name}_${out_normal_label}"

  expected=()
  if [ "$mode" = "all" ] || [ "$mode" = "germline" ]; then
    expected+=(
      "${germline_dir}/${name}_${gdna1_vcf_extension}"
      "${germline_dir}/${name}_${gdna2_vcf_extension}"
      "${germline_dir}/${name}_${gdna3_vcf_extension}"
      "${germline_dir}/${name}_${gdna4_vcf_extension}"
    )
  fi
  if [ "$mode" = "all" ] || [ "$mode" = "rna" ]; then
    rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
    expected+=(
      "${outdir}/${name}_${rna1_vcf_extension}"
      "${outdir}/${name}_${rna2_labeled_vcf_extension}"
      "${outdir}/${name}_${rna3_knownsites_vcf_extension}"
      "${outdir}/${name}_${rna4_summarised_vcf_extension}"
      "${outdir}/${name}_${rna5_qced_vcf_extension}"
      "${outdir}/${name}_${rna6_merged_vcf_extension}"
      "${bamdir}/${name}_${rna_dir_label}/${name}_${rna7_smfixed_bam_suffix}"
      "${outdir}/${name}_${rna7_phased_vcf_extension}"
    )
  fi

  for out in "${expected[@]}"; do
    if [ ! -f "$out" ]; then
      printf "%s\tMISSING\t%s\n" "$name" "$out"
      failed=1
      continue
    fi
    if [ ! -s "$out" ]; then
      printf "%s\tEMPTY\t%s\n" "$name" "$out"
      failed=1
    fi
  done
done < "$samples"

if [ "$checked_samples" -eq 0 ]; then
  echo "ERROR: no matching samples found in $samples" >&2
  exit 1
fi

if [ "$failed" -eq 0 ]; then
  if [ "$mode" = "splicing" ]; then
    echo "All expected splicing step outputs are complete."
  else
    echo "All expected files found and non-empty."
  fi
else
  exit 1
fi
