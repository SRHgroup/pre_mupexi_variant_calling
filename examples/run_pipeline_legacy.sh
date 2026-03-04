#!/usr/bin/bash
set -euo pipefail

# Legacy wrapper for old numeric-step pipeline (4.x / 2.0-3.0).
# Copy to project folder and edit REPO/CONFIG.
REPO="/home/projects/SRHgroup/apps/pre_mupexi_variant_calling_old"
CONFIG="/home/projects/SRHgroup/projects/PemBOv_trial/bin/rnadnavar/CONFIG"

cmd="${1:-help}"
shift || true

if [ "$cmd" = "-h" ] || [ "$cmd" = "--help" ] || [ "$cmd" = "help" ]; then
  cat <<USAGE
Usage:
  $0 rna [PATIENT] [-f]
  $0 germline [PATIENT] [-f]
  $0 all [PATIENT] [-f]
  $0 step <4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0> [PATIENT] [-f]
  $0 check [PATIENT] [all|rna|germline]
  $0 check-step <4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0> [PATIENT]
  $0 watch-step <4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0> [PATIENT] [INTERVAL_SEC]
  $0 sync
  $0 show-config
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
  : "${bamdir:?CONFIG must define bamdir}"
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

legacy_step_expected_output() {
  local patient="$1"
  local step="$2"
  local out_normal out_rna outdir germdir rna_dir_label
  out_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  out_rna="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  outdir="${vcfdir}/${patient}_${out_rna}_vs_${patient}_${out_normal}"
  germdir="${vcfdir}/${patient}_${out_normal}"
  rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"

  case "$step" in
    2.0) printf '%s\n' "${germdir}/${patient}_${output_extension_20}" ;;
    2.0.1) printf '%s\n' "${germdir}/${patient}_${output_extension_201}" ;;
    2.0.2) printf '%s\n' "${germdir}/${patient}_${output_extension_202}" ;;
    3.0) printf '%s\n' "${germdir}/${patient}_${output_extension_30}" ;;
    4.1) printf '%s\n' "${outdir}/${patient}_${rna_only_vcf_extension}" ;;
    4.2) printf '%s\n' "${outdir}/${patient}_${filtered_edit_labeled_vcf_extension}" ;;
    4.3) printf '%s\n' "${outdir}/${patient}_${annot_vcf_extension}" ;;
    4.4) printf '%s\n' "${outdir}/${patient}_${rna_summarised_vcf_extension}" ;;
    4.5) printf '%s\n' "${outdir}/${patient}_${rna_vcf_knownsites_extension}" ;;
    4.6) printf '%s\n' "${outdir}/${patient}_${merged_vcf_extension}" ;;
    4.7.0) printf '%s\n' "${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna_bam_smfixed_suffix}" ;;
    4.7) printf '%s\n' "${outdir}/${patient}_${phased_vcf_extension}" ;;
    *) return 1 ;;
  esac
}

legacy_step_prefix() {
  case "$1" in
    4.1) printf '%s\n' "4.1_OnlyRnaVcf" ;;
    4.2) printf '%s\n' "4.2_FilterEditSignature" ;;
    4.3) printf '%s\n' "4.3_AnnotateKnownSites" ;;
    4.4) printf '%s\n' "4.4_SummariseRnaMetrics" ;;
    4.5) printf '%s\n' "4.5_FilterByAfDpAr" ;;
    4.6) printf '%s\n' "4.6_MergeDnaRnaVcfs" ;;
    4.7.0) printf '%s\n' "4.7.0_FixRnaBamReadGroups" ;;
    4.7) printf '%s\n' "4.7.1_GenotypeAndPhaseMergedVcf" ;;
    2.0) printf '%s\n' "2.0_HaplotypeCaller" ;;
    2.0.1) printf '%s\n' "2.0.1_FilterGermline" ;;
    2.0.2) printf '%s\n' "2.0.2_SelectVariants" ;;
    3.0) printf '%s\n' "3.0_FilterGermlineByAdjacency" ;;
    *) return 1 ;;
  esac
}

check_step_outputs() {
  local step="$1"
  local selected="${2:-}"
  case "$step" in
    4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0) ;;
    *) echo "Unknown step for check-step: $step" >&2; exit 1 ;;
  esac

  load_config
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

    out="$(legacy_step_expected_output "$patient" "$step")"
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

pbs_state_for_jobid() {
  local jobid="$1"
  local state
  state="$(qstat -f "$jobid" 2>/dev/null | awk '/job_state =/{print $3; exit}')"
  case "$state" in
    R|E) printf '%s\n' "RUNNING" ;;
    Q|H|W|T|S) printf '%s\n' "QUEUED" ;;
    C) printf '%s\n' "COMPLETED" ;;
    *) printf '%s\n' "" ;;
  esac
}

watch_step_outputs() {
  local step="$1"
  local selected="${2:-}"
  local interval="${3:-5}"
  case "$step" in
    4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0) ;;
    *) echo "Unknown step for watch-step: $step" >&2; exit 1 ;;
  esac
  [[ "$interval" =~ ^[0-9]+$ ]] || { echo "INTERVAL_SEC must be integer" >&2; exit 1; }
  [ "$interval" -gt 0 ] || { echo "INTERVAL_SEC must be > 0" >&2; exit 1; }

  load_config
  local prefix logdir spinner_idx=0
  prefix="$(legacy_step_prefix "$step")"
  logdir="${vcfdir}/${prefix}.logs_and_reports/logs"
  sp='|/-\\'

  while :; do
    spinner_char="${sp:$spinner_idx:1}"
    spinner_idx=$(( (spinner_idx + 1) % 4 ))

    total=0; done_ok=0; running=0; queued=0; failed=0; not_submitted=0
    rows=()
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
      total=$((total + 1))

      out="$(legacy_step_expected_output "$patient" "$step")"
      marker="${logdir}/submitted.${prefix}.${patient}.jobid"
      status=""
      detail=""

      if [ -f "$out" ] && [ -s "$out" ]; then
        status="DONE"
        detail="$out"
        done_ok=$((done_ok + 1))
      else
        if [ -f "$marker" ]; then
          jobid="$(head -n1 "$marker" 2>/dev/null || true)"
          if [ -n "$jobid" ]; then
            pbs_state="$(pbs_state_for_jobid "$jobid")"
            if [ "$pbs_state" = "RUNNING" ]; then
              status="RUNNING"
              detail="$jobid"
              running=$((running + 1))
            elif [ "$pbs_state" = "QUEUED" ]; then
              status="QUEUED"
              detail="$jobid"
              queued=$((queued + 1))
            else
              status="FAILED_OR_MISSING"
              detail="$out"
              failed=$((failed + 1))
            fi
          else
            status="NOT_SUBMITTED"
            detail="$out"
            not_submitted=$((not_submitted + 1))
          fi
        else
          status="NOT_SUBMITTED"
          detail="$out"
          not_submitted=$((not_submitted + 1))
        fi
      fi

      rows+=("$(printf '%-14s %-18s %s' "$patient" "$status" "$detail")")
    done < "$samples"

    clear
    echo "[$spinner_char] watch-step $step   interval=${interval}s   time=$(date '+%F %T')"
    echo "Totals: total=${total} done=${done_ok} running=${running} queued=${queued} failed=${failed} not_submitted=${not_submitted}"
    echo
    printf '%-14s %-18s %s\n' "PATIENT" "STATUS" "DETAIL"
    printf '%-14s %-18s %s\n' "-------" "------" "------"
    for r in "${rows[@]}"; do
      printf '%s\n' "$r"
    done
    echo
    echo "Ctrl+C to stop."
    sleep "$interval"
  done
}

case "$cmd" in
  rna)       run_make run_all_rna "${1:-}" ;;
  germline)  run_make run_all_germline "${1:-}" ;;
  all)       run_make run_all "${1:-}" ;;
  step)
    step_name="${1:-}"
    sample="${2:-}"
    case "$step_name" in
      4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0)
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
  watch-step)
    step_name="${1:-}"
    sample="${2:-}"
    interval="${3:-5}"
    watch_step_outputs "$step_name" "$sample" "$interval"
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
