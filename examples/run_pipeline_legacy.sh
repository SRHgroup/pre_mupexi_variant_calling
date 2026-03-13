#!/usr/bin/bash
set -euo pipefail

# Legacy-naming wrapper on top of the current pipeline implementation.
# Keeps numeric step names (4.x / 2.0-3.0) but executes targets in the new repo.
# Copy to project folder and edit REPO/CONFIG if needed.
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
  $0 mupexi [PATIENT] [--outdir DIR] [--run-fusions] [--fusion-only] [--hla HLA_STRING] [--expr EXPR_TSV] [--fusion FUSION_ARRIBA_TSV] [--nodes N] [--ppn N] [--mem SIZE] [--walltime HH:MM:SS] [-f] [--skip-running]
  $0 all [PATIENT] [-f]
  $0 research rna-clusters [PATIENT] [--outdir DIR] [--max-distance N] [--min-cluster-size N] [--min-alt-count N] [-f] [--skip-running]
  $0 research samecopy-stats [PATIENT] [--outfile FILE] [--window N] [-f] [--skip-running]
  $0 research variant-table [PATIENT] [--outdir DIR] [-f] [--skip-running]
  $0 research strand-blacklist [PATIENT] [--outdir DIR] [--protocol NAME] [--min-mapq N] [--min-baseq N] [--min-expected-frac X] [-f] [--skip-running]
  $0 research run_mosdepth_overlap [PATIENT] [--outdir DIR] [--depth-threshold N] [--region-bin-size N]
  $0 step <4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0> [PATIENT] [-f]
  $0 check [PATIENT] [all|rna|germline]
  $0 check-step <4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0> [PATIENT]
  $0 watch-step <4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0> [PATIENT] [INTERVAL_SEC]
  $0 sync
  $0 show-config
  (append --skip-running to submit commands to avoid re-submitting active jobs)
USAGE
  exit 0
fi

force=0
skip_running=0
filtered_args=()
while [ $# -gt 0 ]; do
  case "${1:-}" in
    -f|--force) force=1 ;;
    --skip-running) skip_running=1 ;;
    *) filtered_args+=("$1") ;;
  esac
  shift
done
set -- "${filtered_args[@]}"
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

legacy_target() {
  case "$1" in
    4.1) printf '%s\n' "runrna1" ;;
    4.2) printf '%s\n' "runrna2" ;;
    4.3) printf '%s\n' "runrna3" ;;
    4.4) printf '%s\n' "runrna4" ;;
    4.5) printf '%s\n' "runrna5" ;;
    4.6) printf '%s\n' "runrna6" ;;
    4.7.0) printf '%s\n' "runrna7.0" ;;
    4.7) printf '%s\n' "runrna7" ;;
    2.0) printf '%s\n' "rungdna1" ;;
    2.0.1) printf '%s\n' "rungdna2" ;;
    2.0.2) printf '%s\n' "rungdna3" ;;
    3.0) printf '%s\n' "rungdna4" ;;
    *) return 1 ;;
  esac
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

run_research_samecopy_stats() {
  local sample="${1:-}"
  local outfile="${2:-}"
  local window="${3:-33}"

  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_samecopy_stats CONFIG="$CONFIG" SAMPLE="$sample" OUTFILE="$outfile" WINDOW="$window" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_samecopy_stats CONFIG="$CONFIG" OUTFILE="$outfile" WINDOW="$window" $force_arg $skip_running_arg
  fi
}

run_research_variant_table() {
  local sample="${1:-}"
  local outdir="${2:-}"

  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_variant_table CONFIG="$CONFIG" SAMPLE="$sample" OUTDIR="$outdir" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_variant_table CONFIG="$CONFIG" OUTDIR="$outdir" $force_arg $skip_running_arg
  fi
}

run_research_strand_blacklist() {
  local sample="${1:-}"
  local outdir="${2:-}"
  local protocol="${3:-}"
  local min_mapq="${4:-20}"
  local min_baseq="${5:-20}"
  local min_expected_frac="${6:-0.8}"

  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_strand_blacklist CONFIG="$CONFIG" SAMPLE="$sample" OUTDIR="$outdir" PROTOCOL="$protocol" MIN_MAPQ="$min_mapq" MIN_BASEQ="$min_baseq" MIN_EXPECTED_FRAC="$min_expected_frac" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_strand_blacklist CONFIG="$CONFIG" OUTDIR="$outdir" PROTOCOL="$protocol" MIN_MAPQ="$min_mapq" MIN_BASEQ="$min_baseq" MIN_EXPECTED_FRAC="$min_expected_frac" $force_arg $skip_running_arg
  fi
}

run_research_mosdepth_overlap() {
  local sample="${1:-}"
  local outdir="${2:-}"
  local depth_threshold="${3:-10}"
  local region_bin_size="${4:-250000}"

  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_mosdepth_overlap CONFIG="$CONFIG" SAMPLE="$sample" OUTDIR="$outdir" DEPTH_THRESHOLD="$depth_threshold" REGION_BIN_SIZE="$region_bin_size" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_mosdepth_overlap CONFIG="$CONFIG" OUTDIR="$outdir" DEPTH_THRESHOLD="$depth_threshold" REGION_BIN_SIZE="$region_bin_size" $force_arg $skip_running_arg
  fi
}

run_mupexi() {
  local sample="${1:-}"
  local outdir="${2:-}"
  local run_fusions="${3:-0}"
  local fusion_only="${4:-0}"
  local hla="${5:-}"
  local expr="${6:-}"
  local fusion="${7:-}"
  local nodes="${8:-}"
  local ppn="${9:-}"
  local mem="${10:-}"
  local walltime="${11:-}"
  local fusion_arg=""
  local fusion_only_arg=""
  local hla_arg=""
  local expr_arg=""
  local fusion_path_arg=""
  local nodes_arg=""
  local ppn_arg=""
  local mem_arg=""
  local walltime_arg=""
  if [ "$run_fusions" = "1" ]; then
    fusion_arg="RUN_FUSIONS=1"
  fi
  if [ "$fusion_only" = "1" ]; then
    fusion_only_arg="FUSION_ONLY=1"
  fi
  if [ -n "$hla" ]; then hla_arg="HLA=$hla"; fi
  if [ -n "$expr" ]; then expr_arg="EXPR=$expr"; fi
  if [ -n "$fusion" ]; then fusion_path_arg="FUSION=$fusion"; fi
  if [ -n "$nodes" ]; then nodes_arg="MUPEXI_NODES=$nodes"; fi
  if [ -n "$ppn" ]; then ppn_arg="MUPEXI_PPN=$ppn"; fi
  if [ -n "$mem" ]; then mem_arg="MUPEXI_MEM=$mem"; fi
  if [ -n "$walltime" ]; then walltime_arg="MUPEXI_WALLTIME=$walltime"; fi
  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_mupexi CONFIG="$CONFIG" SAMPLE="$sample" OUTDIR="$outdir" $fusion_arg $fusion_only_arg $hla_arg $expr_arg $fusion_path_arg $nodes_arg $ppn_arg $mem_arg $walltime_arg $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_mupexi CONFIG="$CONFIG" OUTDIR="$outdir" $fusion_arg $fusion_only_arg $hla_arg $expr_arg $fusion_path_arg $nodes_arg $ppn_arg $mem_arg $walltime_arg $force_arg $skip_running_arg
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

legacy_step_input_status() {
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
    2.0)
      local dna_bam="${bamdir}/${patient}_${out_normal}/${patient}_${dna_bam_suffix:-md.bam}"
      if [ -f "$dna_bam" ] && [ -s "$dna_bam" ]; then
        printf 'INPUT_OK\t%s\n' "$dna_bam"
      else
        printf 'NO_INPUT\t%s\n' "$dna_bam"
      fi
      ;;
    2.0.1)
      local in_201="${germdir}/${patient}_${output_extension_20}"
      if [ -f "$in_201" ] && [ -s "$in_201" ]; then
        printf 'INPUT_OK\t%s\n' "$in_201"
      else
        printf 'NO_INPUT\t%s\n' "$in_201"
      fi
      ;;
    2.0.2)
      local in_202="${germdir}/${patient}_${output_extension_201}"
      if [ -f "$in_202" ] && [ -s "$in_202" ]; then
        printf 'INPUT_OK\t%s\n' "$in_202"
      else
        printf 'NO_INPUT\t%s\n' "$in_202"
      fi
      ;;
    3.0)
      local in_30="${germdir}/${patient}_${output_extension_202}"
      if [ -f "$in_30" ] && [ -s "$in_30" ]; then
        printf 'INPUT_OK\t%s\n' "$in_30"
      else
        printf 'NO_INPUT\t%s\n' "$in_30"
      fi
      ;;
    4.1)
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
    4.2)
      local in_42="${outdir}/${patient}_${rna_only_vcf_extension}"
      if [ -f "$in_42" ] && [ -s "$in_42" ]; then
        printf 'INPUT_OK\t%s\n' "$in_42"
      else
        printf 'NO_INPUT\t%s\n' "$in_42"
      fi
      ;;
    4.3)
      local in_43="${outdir}/${patient}_${filtered_edit_labeled_vcf_extension}"
      if [ -f "$in_43" ] && [ -s "$in_43" ]; then
        printf 'INPUT_OK\t%s\n' "$in_43"
      else
        printf 'NO_INPUT\t%s\n' "$in_43"
      fi
      ;;
    4.4)
      local in_44="${outdir}/${patient}_${annot_vcf_extension}"
      if [ -f "$in_44" ] && [ -s "$in_44" ]; then
        printf 'INPUT_OK\t%s\n' "$in_44"
      else
        printf 'NO_INPUT\t%s\n' "$in_44"
      fi
      ;;
    4.5)
      local in_45="${outdir}/${patient}_${rna_summarised_vcf_extension}"
      if [ -f "$in_45" ] && [ -s "$in_45" ]; then
        printf 'INPUT_OK\t%s\n' "$in_45"
      else
        printf 'NO_INPUT\t%s\n' "$in_45"
      fi
      ;;
    4.6)
      local in_46="${outdir}/${patient}_${rna_vcf_knownsites_extension}"
      if [ -f "$in_46" ] && [ -s "$in_46" ]; then
        printf 'INPUT_OK\t%s\n' "$in_46"
      else
        printf 'NO_INPUT\t%s\n' "$in_46"
      fi
      ;;
    4.7.0)
      local in_470="${bamdir}/${patient}_${rna_dir_label}/${patient}_${dna_bam_suffix:-md.bam}"
      if [ -f "$in_470" ] && [ -s "$in_470" ]; then
        printf 'INPUT_OK\t%s\n' "$in_470"
      else
        printf 'NO_INPUT\t%s\n' "$in_470"
      fi
      ;;
    4.7)
      local in_47_1="${outdir}/${patient}_${merged_vcf_extension}"
      local in_47_2="${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna_bam_smfixed_suffix}"
      if [ -f "$in_47_1" ] && [ -s "$in_47_1" ] && [ -f "$in_47_2" ] && [ -s "$in_47_2" ]; then
        printf 'INPUT_OK\tVCF:%s ; BAM:%s\n' "$in_47_1" "$in_47_2"
      elif [ ! -f "$in_47_1" ] || [ ! -s "$in_47_1" ]; then
        printf 'NO_INPUT\tVCF:%s\n' "$in_47_1"
      else
        printf 'NO_INPUT\tBAM:%s\n' "$in_47_2"
      fi
      ;;
    *)
      printf 'NO_INPUT\tunknown-step\n'
      ;;
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
      input_info="$(legacy_step_input_status "$patient" "$step")"
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
  dna-only)  run_make run_dna_only "${1:-}" ;;
  mupexi)
    sample=""
    outdir=""
    run_fusions="0"
    fusion_only="0"
    hla=""
    expr=""
    fusion=""
    mupexi_nodes=""
    mupexi_ppn=""
    mupexi_mem=""
    mupexi_walltime=""
    if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
      sample="$1"
      shift
    fi
    while [ $# -gt 0 ]; do
      case "${1:-}" in
        --outdir|-o) outdir="${2:-}"; shift 2 ;;
        --run-fusions) run_fusions="1"; shift ;;
        --fusion-only) fusion_only="1"; run_fusions="1"; shift ;;
        --hla) hla="${2:-}"; shift 2 ;;
        --expr) expr="${2:-}"; shift 2 ;;
        --fusion) fusion="${2:-}"; shift 2 ;;
        --nodes) mupexi_nodes="${2:-}"; shift 2 ;;
        --ppn) mupexi_ppn="${2:-}"; shift 2 ;;
        --mem) mupexi_mem="${2:-}"; shift 2 ;;
        --walltime) mupexi_walltime="${2:-}"; shift 2 ;;
        *) echo "Unknown mupexi option: $1" >&2; exit 1 ;;
      esac
    done
    run_mupexi "$sample" "$outdir" "$run_fusions" "$fusion_only" "$hla" "$expr" "$fusion" "$mupexi_nodes" "$mupexi_ppn" "$mupexi_mem" "$mupexi_walltime"
    ;;
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
      samecopy-stats)
        sample=""
        outfile=""
        window="33"
        if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
          sample="$1"
          shift
        fi
        while [ $# -gt 0 ]; do
          case "${1:-}" in
            --outfile|-o) outfile="${2:-}"; shift 2 ;;
            --window) window="${2:-}"; shift 2 ;;
            *) echo "Unknown research option: $1" >&2; exit 1 ;;
          esac
        done
        run_research_samecopy_stats "$sample" "$outfile" "$window"
        ;;
      variant-table)
        sample=""
        outdir=""
        if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
          sample="$1"
          shift
        fi
        while [ $# -gt 0 ]; do
          case "${1:-}" in
            --outdir|-o) outdir="${2:-}"; shift 2 ;;
            *) echo "Unknown research option: $1" >&2; exit 1 ;;
          esac
        done
        run_research_variant_table "$sample" "$outdir"
        ;;
      strand-blacklist)
        sample=""
        outdir=""
        protocol=""
        min_mapq="20"
        min_baseq="20"
        min_expected_frac="0.8"
        if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
          sample="$1"
          shift
        fi
        while [ $# -gt 0 ]; do
          case "${1:-}" in
            --outdir|-o) outdir="${2:-}"; shift 2 ;;
            --protocol) protocol="${2:-}"; shift 2 ;;
            --min-mapq) min_mapq="${2:-}"; shift 2 ;;
            --min-baseq) min_baseq="${2:-}"; shift 2 ;;
            --min-expected-frac) min_expected_frac="${2:-}"; shift 2 ;;
            *) echo "Unknown research option: $1" >&2; exit 1 ;;
          esac
        done
        run_research_strand_blacklist "$sample" "$outdir" "$protocol" "$min_mapq" "$min_baseq" "$min_expected_frac"
        ;;
      run_mosdepth_overlap|mosdepth-overlap)
        sample=""
        outdir=""
        depth_threshold="10"
        region_bin_size="250000"
        if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
          sample="$1"
          shift
        fi
        while [ $# -gt 0 ]; do
          case "${1:-}" in
            --outdir|-o) outdir="${2:-}"; shift 2 ;;
            --depth-threshold) depth_threshold="${2:-}"; shift 2 ;;
            --region-bin-size) region_bin_size="${2:-}"; shift 2 ;;
            *) echo "Unknown research option: $1" >&2; exit 1 ;;
          esac
        done
        run_research_mosdepth_overlap "$sample" "$outdir" "$depth_threshold" "$region_bin_size"
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
      4.1|4.2|4.3|4.4|4.5|4.6|4.7.0|4.7|2.0|2.0.1|2.0.2|3.0)
        target="$(legacy_target "$step_name")"
        run_make "$target" "$sample"
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
