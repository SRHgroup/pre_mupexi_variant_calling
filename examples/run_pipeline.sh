#!/usr/bin/bash
set -euo pipefail

# Copy this file to your project folder, then edit these two paths:
REPO="/home/projects/SRHgroup/apps/pre_mupexi_variant_calling"
CONFIG="/home/projects/SRHgroup/projects/PemBOv_trial/bin/rnadnavar/CONFIG"

cmd="${1:-help}"
shift || true

if [ "$cmd" = "-h" ] || [ "$cmd" = "--help" ] || [ "$cmd" = "help" ]; then
  cat <<USAGE
Usage:
  $0 rna [PATIENT] [-f]
  $0 germline [PATIENT] [-f]
  $0 all [PATIENT] [-f]
  $0 step <rna1|rna2|rna3|rna4|rna5|rna6|rna7|gdna1|gdna2|gdna3|gdna4> [PATIENT] [-f]
  $0 check [PATIENT] [all|rna|germline]
  $0 check-step <rna1|rna2|rna3|rna4|rna5|rna6|rna7|gdna1|gdna2|gdna3|gdna4> [PATIENT]
  $0 sync
  $0 show-config

Examples:
  $0 show-config
  $0 sync
  $0 rna
  $0 rna 01-CH-L
  $0 germline 01-CH-L
  $0 all 01-CH-L
  $0 step rna4 01-CH-L -f
  $0 step gdna2 01-CH-L
  $0 check
  $0 check 01-CH-L rna
  $0 check-step rna5
  $0 check-step gdna4 01-CH-L
USAGE
  exit 0
fi

force=0
if [[ " ${*:-} " == *" -f "* ]] || [[ " ${*:-} " == *" --force "* ]]; then
  force=1
  set -- "${@/-f/}"
  set -- "${@/--force/}"
fi
force_arg=""
if [ "$force" -eq 1 ]; then
  force_arg="FORCE=1"
fi

run_make() {
  local target="$1"
  local sample="${2:-}"
  if [ -n "$sample" ]; then
    make -C "$REPO" "$target" CONFIG="$CONFIG" SAMPLE="$sample" $force_arg
  else
    make -C "$REPO" "$target" CONFIG="$CONFIG" $force_arg
  fi
}

load_config() {
  [ -f "$CONFIG" ] || { echo "ERROR: CONFIG not found: $CONFIG" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$CONFIG"
  : "${samples:?CONFIG must define samples}"
  : "${vcfdir:?CONFIG must define vcfdir}"
}

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label:-DNA_NORMAL}"
    "${dna_tumor_label:-DNA_TUMOR}"
    "${rna_tumor_label:-RNA_TUMOR}"
    "${out_dna_normal_label:-DNA_NORMAL}"
    "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
    "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do
    value="${value%_${label}}"
  done
  printf '%s\n' "$value"
}

step_expected_output() {
  local patient="$1"
  local step="$2"
  local out_normal out_rna outdir germdir
  out_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  out_rna="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  outdir="${vcfdir}/${patient}_${out_rna}_vs_${patient}_${out_normal}"
  germdir="${vcfdir}/${patient}_${out_normal}"
  case "$step" in
    gdna1) printf '%s\n' "${germdir}/${patient}_${gdna1_vcf_extension}" ;;
    gdna2) printf '%s\n' "${germdir}/${patient}_${gdna2_vcf_extension}" ;;
    gdna3) printf '%s\n' "${germdir}/${patient}_${gdna3_vcf_extension}" ;;
    gdna4) printf '%s\n' "${germdir}/${patient}_${gdna4_vcf_extension}" ;;
    rna1) printf '%s\n' "${outdir}/${patient}_${rna1_vcf_extension}" ;;
    rna2) printf '%s\n' "${outdir}/${patient}_${rna2_labeled_vcf_extension}" ;;
    rna3) printf '%s\n' "${outdir}/${patient}_${rna3_knownsites_vcf_extension}" ;;
    rna4) printf '%s\n' "${outdir}/${patient}_${rna4_summarised_vcf_extension}" ;;
    rna5) printf '%s\n' "${outdir}/${patient}_${rna5_qced_vcf_extension}" ;;
    rna6) printf '%s\n' "${outdir}/${patient}_${rna6_merged_vcf_extension}" ;;
    rna7) printf '%s\n' "${outdir}/${patient}_${rna7_phased_vcf_extension}" ;;
    *) return 1 ;;
  esac
}

check_step_outputs() {
  local step="$1"
  local selected="${2:-}"
  case "$step" in
    rna1|rna2|rna3|rna4|rna5|rna6|rna7|gdna1|gdna2|gdna3|gdna4) ;;
    *) echo "Unknown step for check-step: $step" >&2; exit 1 ;;
  esac

  load_config
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
  : "${rna7_phased_vcf_extension:?CONFIG must define rna7_phased_vcf_extension}"

  local failed=0
  declare -A seen=()
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    case "$line" in [[:space:]]*'#'*) continue ;; esac
    sample_id=$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')
    patient=$(sample_base_name "$sample_id")
    [ -n "$patient" ] || continue
    [ -n "${seen[$patient]:-}" ] && continue
    seen["$patient"]=1
    if [ -n "$selected" ] && [ "$selected" != "$patient" ] && [ "$selected" != "$sample_id" ]; then
      continue
    fi

    out="$(step_expected_output "$patient" "$step")"
    if [ -f "$out" ] && [ -s "$out" ]; then
      printf "%s\tYES\t%s\n" "$patient" "$out"
    elif [ -f "$out" ]; then
      printf "%s\tNO\tEMPTY\t%s\n" "$patient" "$out"
      failed=1
    else
      printf "%s\tNO\tMISSING\t%s\n" "$patient" "$out"
      failed=1
    fi
  done < "$samples"

  [ "$failed" -eq 0 ]
}

case "$cmd" in
  rna)       run_make run_all_rna "${1:-}" ;;
  germline)  run_make run_all_germline "${1:-}" ;;
  all)       run_make run_all "${1:-}" ;;
  step)
    step_name="${1:-}"
    sample="${2:-}"
    case "$step_name" in
      rna1|rna2|rna3|rna4|rna5|rna6|rna7|gdna1|gdna2|gdna3|gdna4)
        run_make "run${step_name}" "$sample"
        ;;
      *)
        echo "Unknown step: $step_name" >&2
        exit 1
        ;;
    esac
    ;;
  check)
    sample="${1:-}"
    mode="${2:-all}"
    if [ -n "$sample" ]; then
      make -C "$REPO" check_outputs CONFIG="$CONFIG" SAMPLE="$sample" MODE="$mode"
    else
      make -C "$REPO" check_outputs CONFIG="$CONFIG" MODE="$mode"
    fi
    ;;
  check-step)
    step_name="${1:-}"
    sample="${2:-}"
    check_step_outputs "$step_name" "$sample"
    ;;
  sync)
    git -C "$REPO" pull --ff-only origin main
    ;;
  show-config)
    echo "REPO=$REPO"
    echo "CONFIG=$CONFIG"
    ;;
  *)
    echo "Unknown command: $cmd" >&2
    echo "Use: $0 --help" >&2
    exit 1
    ;;
esac
