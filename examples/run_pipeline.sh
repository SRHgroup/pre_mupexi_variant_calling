#!/usr/bin/bash
set -euo pipefail

# Copy this file to your project folder, then edit these two paths:
REPO="/home/projects/SRHgroup/apps/pre_mupexi_variant_calling"
CONFIG="/home/projects/SRHgroup/projects/PemBOv_trial/bin/rnadnavar/CONFIG"
PIPELINE_DEFAULTS="${PIPELINE_DEFAULTS:-$REPO/pipeline_defaults/toolchain.defaults.sh}"

cmd="${1:-help}"
shift || true

if [ "$cmd" = "-h" ] || [ "$cmd" = "--help" ] || [ "$cmd" = "help" ]; then
  cat <<USAGE
Usage:
  $0 rna [PATIENT] [-f]
  $0 germline [PATIENT] [-f]
  $0 dna-only [PATIENT] [-f]
  $0 all [PATIENT] [-f]
  $0 research rna-clusters [PATIENT] [--outdir DIR] [--max-distance N] [--min-cluster-size N] [--min-alt-count N] [-f] [--skip-running]
  $0 step <rna1|rna2|rna3|rna4|rna5|rna6|rna7.0|rna7|gdna1|gdna2|gdna3|gdna4> [PATIENT] [-f]
  $0 check [PATIENT] [all|rna|germline]
  $0 check-step <rna1|rna2|rna3|rna4|rna5|rna6|rna7.0|rna7|gdna1|gdna2|gdna3|gdna4> [PATIENT]
  $0 watch-step <rna1|rna2|rna3|rna4|rna5|rna6|rna7.0|rna7|gdna1|gdna2|gdna3|gdna4> [PATIENT] [INTERVAL_SEC]
  $0 sync
  $0 show-config
  (append --skip-running to submit commands to avoid re-submitting active jobs)

Examples:
  $0 show-config
  $0 sync
  $0 rna
  $0 rna 01-CH-L
  $0 germline 01-CH-L
  $0 dna-only 01-CH-L
  $0 all 01-CH-L
  $0 step rna4 01-CH-L -f
  $0 step rna7 01-CH-L --skip-running
  $0 step rna7.0 01-CH-L -f
  $0 step rna7 01-CH-L -f
  $0 step gdna2 01-CH-L
  $0 research rna-clusters --outdir /home/projects/SRHgroup/projects/SingelCell_Bladder/data/rna/rnadnavar/editing_clusters33
  $0 research rna-clusters Pat96 --max-distance 33 --min-alt-count 10
  $0 check
  $0 check 01-CH-L rna
  $0 check-step rna5
  $0 check-step gdna4 01-CH-L
  $0 watch-step rna7
  $0 watch-step gdna3 01-CH-L 10
USAGE
  exit 0
fi

force=0
skip_running=0
if [[ " ${*:-} " == *" -f "* ]] || [[ " ${*:-} " == *" --force "* ]]; then
  force=1
  set -- "${@/-f/}"
  set -- "${@/--force/}"
fi
if [[ " ${*:-} " == *" --skip-running "* ]]; then
  skip_running=1
  set -- "${@/--skip-running/}"
fi
force_arg=""
if [ "$force" -eq 1 ]; then
  force_arg="FORCE=1"
fi
skip_running_arg=""
if [ "$skip_running" -eq 1 ]; then
  skip_running_arg="SKIP_RUNNING=1"
fi

run_make() {
  local target="$1"
  local sample="${2:-}"
  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" "$target" CONFIG="$CONFIG" SAMPLE="$sample" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" "$target" CONFIG="$CONFIG" $force_arg $skip_running_arg
  fi
}

run_research_rna_clusters() {
  local sample="${1:-}"
  local outdir="${2:-}"
  local max_distance="${3:-33}"
  local min_cluster_size="${4:-2}"
  local min_alt_count="${5:-10}"

  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_rna_clusters CONFIG="$CONFIG" SAMPLE="$sample" OUTDIR="$outdir" MAX_DISTANCE="$max_distance" MIN_CLUSTER_SIZE="$min_cluster_size" MIN_ALT_COUNT="$min_alt_count" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_rna_clusters CONFIG="$CONFIG" OUTDIR="$outdir" MAX_DISTANCE="$max_distance" MIN_CLUSTER_SIZE="$min_cluster_size" MIN_ALT_COUNT="$min_alt_count" $force_arg $skip_running_arg
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

resolve_patient_placeholder() {
  local template="$1"
  local patient="$2"
  local token='{patient}'
  printf '%s\n' "${template//"$token"/$patient}"
}

pick_first_existing() {
  local path
  for path in "$@"; do
    [ -n "$path" ] || continue
    if [ -f "$path" ]; then
      printf '%s\n' "$path"
      return 0
    fi
  done
  return 1
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
    rna7.0)
      rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
      printf '%s\n' "${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna7_smfixed_bam_suffix}"
      ;;
    rna7) printf '%s\n' "${outdir}/${patient}_${rna7_phased_vcf_extension}" ;;
    *) return 1 ;;
  esac
}

step_input_status() {
  local patient="$1"
  local step="$2"
  local out_normal out_rna dna_label rna_dir_label outdir germdir
  local rna_vcf_dir dna_vcf_dir srna_pref sdna_pref srna sdna
  out_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  out_rna="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
  rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
  outdir="${vcfdir}/${patient}_${out_rna}_vs_${patient}_${out_normal}"
  germdir="${vcfdir}/${patient}_${out_normal}"

  case "$step" in
    gdna1)
      local dna_bam="${bamdir}/${patient}_${out_normal}/${patient}_${dna_bam_suffix:-md.bam}"
      if [ -f "$dna_bam" ] && [ -s "$dna_bam" ]; then
        printf 'INPUT_OK\t%s\n' "$dna_bam"
      else
        printf 'NO_INPUT\t%s\n' "$dna_bam"
      fi
      ;;
    gdna2)
      local in_gdna2="${germdir}/${patient}_${gdna1_vcf_extension}"
      if [ -f "$in_gdna2" ] && [ -s "$in_gdna2" ]; then
        printf 'INPUT_OK\t%s\n' "$in_gdna2"
      else
        printf 'NO_INPUT\t%s\n' "$in_gdna2"
      fi
      ;;
    gdna3)
      local in_gdna3="${germdir}/${patient}_${gdna2_vcf_extension}"
      if [ -f "$in_gdna3" ] && [ -s "$in_gdna3" ]; then
        printf 'INPUT_OK\t%s\n' "$in_gdna3"
      else
        printf 'NO_INPUT\t%s\n' "$in_gdna3"
      fi
      ;;
    gdna4)
      local in_gdna4="${germdir}/${patient}_${gdna3_vcf_extension}"
      if [ -f "$in_gdna4" ] && [ -s "$in_gdna4" ]; then
        printf 'INPUT_OK\t%s\n' "$in_gdna4"
      else
        printf 'NO_INPUT\t%s\n' "$in_gdna4"
      fi
      ;;
    rna1)
      source_rna_mutect2_vcf_extension="${source_rna_mutect2_vcf_extension:-${out_rna}_vs_{patient}_${out_normal}.mutect2.filtered.vcf.gz}"
      source_dna_mutect2_vcf_extension="${source_dna_mutect2_vcf_extension:-${dna_label}_vs_{patient}_${out_normal}.mutect2.filtered.vcf.gz}"
      source_rna_mutect2_vcf_extension="$(resolve_patient_placeholder "$source_rna_mutect2_vcf_extension" "$patient")"
      source_dna_mutect2_vcf_extension="$(resolve_patient_placeholder "$source_dna_mutect2_vcf_extension" "$patient")"
      source_rna_mutect2_vcf_extension="${source_rna_mutect2_vcf_extension%\}}"
      source_dna_mutect2_vcf_extension="${source_dna_mutect2_vcf_extension%\}}"
      rna_vcf_dir="${vcfdir}/${patient}_${out_rna}_vs_${patient}_${out_normal}"
      dna_vcf_dir="${vcfdir}/${patient}_${dna_label}_vs_${patient}_${out_normal}"
      srna_pref="${rna_vcf_dir}/${patient}_${source_rna_mutect2_vcf_extension}"
      sdna_pref="${dna_vcf_dir}/${patient}_${source_dna_mutect2_vcf_extension}"
      srna="$(pick_first_existing \
        "$srna_pref" \
        "${rna_vcf_dir}/${patient}_${out_rna}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
        "${rna_vcf_dir}/${patient}_${rna_tumor_label:-RNA_TUMOR}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
        "${rna_vcf_dir}/${patient}_RNA_TUMOR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
        "${rna_vcf_dir}/${patient}_RNA_TUMOUR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      )" || srna=""
      sdna="$(pick_first_existing \
        "$sdna_pref" \
        "${dna_vcf_dir}/${patient}_${dna_label}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
        "${dna_vcf_dir}/${patient}_${dna_tumor_label:-DNA_TUMOR}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
        "${dna_vcf_dir}/${patient}_DNA_TUMOR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
        "${dna_vcf_dir}/${patient}_DNA_TUMOUR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      )" || sdna=""
      if [ -z "$srna" ]; then
        printf 'NO_INPUT\tRNA:%s\n' "$srna_pref"
      elif [ -z "$sdna" ]; then
        printf 'NO_INPUT\tDNA:%s\n' "$sdna_pref"
      else
        printf 'INPUT_OK\tRNA:%s ; DNA:%s\n' "$srna" "$sdna"
      fi
      ;;
    rna2)
      local in_rna2="${outdir}/${patient}_${rna1_vcf_extension}"
      if [ -f "$in_rna2" ] && [ -s "$in_rna2" ]; then
        printf 'INPUT_OK\t%s\n' "$in_rna2"
      else
        printf 'NO_INPUT\t%s\n' "$in_rna2"
      fi
      ;;
    rna3)
      local in_rna3="${outdir}/${patient}_${rna2_labeled_vcf_extension}"
      if [ -f "$in_rna3" ] && [ -s "$in_rna3" ]; then
        printf 'INPUT_OK\t%s\n' "$in_rna3"
      else
        printf 'NO_INPUT\t%s\n' "$in_rna3"
      fi
      ;;
    rna4)
      local in_rna4="${outdir}/${patient}_${rna3_knownsites_vcf_extension}"
      if [ -f "$in_rna4" ] && [ -s "$in_rna4" ]; then
        printf 'INPUT_OK\t%s\n' "$in_rna4"
      else
        printf 'NO_INPUT\t%s\n' "$in_rna4"
      fi
      ;;
    rna5)
      local in_rna5="${outdir}/${patient}_${rna4_summarised_vcf_extension}"
      if [ -f "$in_rna5" ] && [ -s "$in_rna5" ]; then
        printf 'INPUT_OK\t%s\n' "$in_rna5"
      else
        printf 'NO_INPUT\t%s\n' "$in_rna5"
      fi
      ;;
    rna6)
      local in_rna6="${outdir}/${patient}_${rna5_qced_vcf_extension}"
      if [ -f "$in_rna6" ] && [ -s "$in_rna6" ]; then
        printf 'INPUT_OK\t%s\n' "$in_rna6"
      else
        printf 'NO_INPUT\t%s\n' "$in_rna6"
      fi
      ;;
    rna7.0)
      local in_rna70="${bamdir}/${patient}_${rna_dir_label}/${patient}_${dna_bam_suffix:-md.bam}"
      if [ -f "$in_rna70" ] && [ -s "$in_rna70" ]; then
        printf 'INPUT_OK\t%s\n' "$in_rna70"
      else
        printf 'NO_INPUT\t%s\n' "$in_rna70"
      fi
      ;;
    rna7)
      local in_rna7_1="${outdir}/${patient}_${rna6_merged_vcf_extension}"
      local in_rna7_2="${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna7_smfixed_bam_suffix}"
      if [ -f "$in_rna7_1" ] && [ -s "$in_rna7_1" ] && [ -f "$in_rna7_2" ] && [ -s "$in_rna7_2" ]; then
        printf 'INPUT_OK\tVCF:%s ; BAM:%s\n' "$in_rna7_1" "$in_rna7_2"
      elif [ ! -f "$in_rna7_1" ] || [ ! -s "$in_rna7_1" ]; then
        printf 'NO_INPUT\tVCF:%s\n' "$in_rna7_1"
      else
        printf 'NO_INPUT\tBAM:%s\n' "$in_rna7_2"
      fi
      ;;
    *)
      printf 'NO_INPUT\tunknown-step\n'
      ;;
  esac
}

step_prefix() {
  case "$1" in
    rna1) printf '%s\n' "rna1_OnlyRnaVcf" ;;
    rna2) printf '%s\n' "rna2_FilterEditSignature" ;;
    rna3) printf '%s\n' "rna3_AnnotateKnownSites" ;;
    rna4) printf '%s\n' "rna4_SummariseRnaMetrics" ;;
    rna5) printf '%s\n' "rna5_FilterByAfDpAr" ;;
    rna6) printf '%s\n' "rna6_MergeDnaRnaVcfs" ;;
    rna7.0) printf '%s\n' "rna7.0_FixRnaBamReadGroups" ;;
    rna7) printf '%s\n' "rna7_GenotypeAndPhaseMergedVcf" ;;
    gdna1) printf '%s\n' "gdna1_HaplotypeCaller" ;;
    gdna2) printf '%s\n' "gdna2_FilterGermline" ;;
    gdna3) printf '%s\n' "gdna3_SelectVariants" ;;
    gdna4) printf '%s\n' "gdna4_FilterGermlineByAdjacency" ;;
    *) return 1 ;;
  esac
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
    rna1|rna2|rna3|rna4|rna5|rna6|rna7.0|rna7|gdna1|gdna2|gdna3|gdna4) ;;
    *) echo "Unknown step for watch-step: $step" >&2; exit 1 ;;
  esac
  [[ "$interval" =~ ^[0-9]+$ ]] || { echo "INTERVAL_SEC must be integer" >&2; exit 1; }
  [ "$interval" -gt 0 ] || { echo "INTERVAL_SEC must be > 0" >&2; exit 1; }

  load_config
  local prefix logdir spinner_idx=0
  prefix="$(step_prefix "$step")"
  logdir="${vcfdir}/${prefix}.logs_and_reports/logs"
  sp='|/-\'

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

      out="$(step_expected_output "$patient" "$step")"
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

check_step_outputs() {
  local step="$1"
  local selected="${2:-}"
  case "$step" in
    rna1|rna2|rna3|rna4|rna5|rna6|rna7.0|rna7|gdna1|gdna2|gdna3|gdna4) ;;
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
  : "${rna7_smfixed_bam_suffix:?CONFIG must define rna7_smfixed_bam_suffix}"
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
      input_info="$(step_input_status "$patient" "$step")"
      input_state="$(printf '%s\n' "$input_info" | awk -F'\t' 'NR==1{print $1}')"
      input_path="$(printf '%s\n' "$input_info" | awk -F'\t' 'NR==1{print $2}')"
      if [ "$input_state" = "INPUT_OK" ]; then
        printf "%s\tNO\tMISSING_OUTPUT\t%s\tINPUT_OK\t%s\n" "$patient" "$out" "$input_path"
      else
        printf "%s\tNO\tMISSING_OUTPUT\t%s\tNO_INPUT\t%s\n" "$patient" "$out" "$input_path"
      fi
      failed=1
    fi
  done < "$samples"

  [ "$failed" -eq 0 ]
}

case "$cmd" in
  rna)       run_make run_all_rna "${1:-}" ;;
  germline)  run_make run_all_germline "${1:-}" ;;
  dna-only)  run_make run_dna_only "${1:-}" ;;
  all)       run_make run_all "${1:-}" ;;
  research)
    task="${1:-}"
    shift || true
    case "$task" in
      rna-clusters)
        sample=""
        outdir=""
        max_distance="33"
        min_cluster_size="2"
        min_alt_count="10"
        if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
          sample="$1"
          shift
        fi
        while [ $# -gt 0 ]; do
          case "${1:-}" in
            --outdir|-o) outdir="${2:-}"; shift 2 ;;
            --max-distance) max_distance="${2:-}"; shift 2 ;;
            --min-cluster-size) min_cluster_size="${2:-}"; shift 2 ;;
            --min-alt-count) min_alt_count="${2:-}"; shift 2 ;;
            *) echo "Unknown research option: $1" >&2; exit 1 ;;
          esac
        done
        run_research_rna_clusters "$sample" "$outdir" "$max_distance" "$min_cluster_size" "$min_alt_count"
        ;;
      *)
        echo "Unknown research task: $task" >&2
        exit 1
        ;;
    esac
    ;;
  step)
    step_name="${1:-}"
    sample="${2:-}"
    case "$step_name" in
      rna1|rna2|rna3|rna4|rna5|rna6|rna7.0|rna7|gdna1|gdna2|gdna3|gdna4)
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
