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
  $0 splicing merge-star-sj [PATIENT] [--root STAR_DIR] [-f] [--dry-run]
  $0 mupexi [PATIENT] [--outdir DIR] [--run-fusions] [--fusion-only] [--hla HLA_STRING] [--expr EXPR_TSV] [--fusion FUSION_ARRIBA_TSV] [--nodes N] [--ppn N] [--mem SIZE] [--walltime HH:MM:SS] [-f] [--skip-running]
  $0 cleanup [PATIENT] [--execute] [--threads N] [--nodes N] [--ppn N] [--mem SIZE] [--walltime HH:MM:SS] [-f] [--skip-running]
  $0 all [PATIENT] [-f]
  $0 research rna-clusters [PATIENT] [--outdir DIR] [--max-distance N] [--min-cluster-size N] [--min-alt-count N] [-f] [--skip-running]
  $0 research samecopy-stats [PATIENT] [--outfile FILE] [--window N] [-f] [--skip-running]
  $0 research variant-table [PATIENT] [--outdir DIR] [-f] [--skip-running]
  $0 research strand-blacklist [PATIENT] [--outdir DIR] [--protocol NAME] [--min-mapq N] [--min-baseq N] [--min-expected-frac X] [-f] [--skip-running]
  $0 research run_mosdepth_overlap [PATIENT] [--outdir DIR] [--depth-threshold N] [--region-bin-size N]
  $0 step <rna1|rna2|rna3|rna4|rna5|rna5.1|rna6|rna7.0|rna7|rna7.1|gdna1|gdna2|gdna3|gdna4> [PATIENT] [-f]
  $0 check [PATIENT] [all|rna|germline]
  $0 check-step <rna1|rna2|rna3|rna4|rna5|rna5.1|rna6|rna7.0|rna7|rna7.1|gdna1|gdna2|gdna3|gdna4> [PATIENT]
  $0 watch-step <rna1|rna2|rna3|rna4|rna5|rna5.1|rna6|rna7.0|rna7|rna7.1|gdna1|gdna2|gdna3|gdna4> [PATIENT] [INTERVAL_SEC]
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
  $0 splicing merge-star-sj Pat21
  $0 splicing merge-star-sj Pat21 --root /path/to/reports/star --dry-run
  $0 mupexi 01-CH-L
  $0 mupexi 01-CH-L --run-fusions
  $0 mupexi 01-CH-L --fusion-only --run-fusions --outdir /path/to/mupexi2_fusions_only
  $0 mupexi 01-CH-L --ppn 8 --mem 48gb --walltime 24:00:00
  $0 mupexi Pat11 --hla HLA-A26:01,HLA-A24:02,HLA-B27:02,HLA-B15:09,HLA-C02:02,HLA-C07:04 --expr /path/to/Patient_11.expr.tsv
  $0 cleanup Pat101
  $0 cleanup Pat101 --execute --ppn 8 --mem 32gb --walltime 24:00:00
  $0 cleanup --execute --skip-running
  $0 all 01-CH-L
  $0 step rna4 01-CH-L -f
  $0 step rna5.1 01-CH-L -f
  $0 step rna7 01-CH-L --skip-running
  $0 step rna7.1 01-CH-L -f
  $0 step rna7.0 01-CH-L -f
  $0 step rna7 01-CH-L -f
  $0 step gdna2 01-CH-L
  $0 research rna-clusters --outdir /home/projects/SRHgroup/projects/SingelCell_Bladder/data/rna/rnadnavar/editing_clusters33
  $0 research rna-clusters Pat96 --max-distance 33 --min-alt-count 10
  $0 research samecopy-stats --outfile /home/projects/SRHgroup/projects/SingelCell_Bladder/data/rna/rnadnavar/samecopy_stats_with_rna.tsv
  $0 research variant-table
  $0 research variant-table Pat11 --outdir /home/projects/SRHgroup/projects/SingelCell_Bladder/data/rna/rnadnavar/variant_tables
  $0 research strand-blacklist
  $0 research strand-blacklist Pat11 --min-expected-frac 0.9
  $0 research run_mosdepth_overlap --depth-threshold 10 --region-bin-size 250000
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

run_research_vep_dedup() {
  local sample="${1:-}"
  local outdir="${2:-}"

  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_vep_dedup CONFIG="$CONFIG" SAMPLE="$sample" OUTDIR="$outdir" $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_research_vep_dedup CONFIG="$CONFIG" OUTDIR="$outdir" $force_arg $skip_running_arg
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

run_splicing_merge_star_sj() {
  local sample="${1:-}"
  local star_root="${2:-}"
  local dry_run="${3:-0}"
  local star_root_arg=""
  local dry_run_arg=""
  if [ -n "$star_root" ]; then star_root_arg="STAR_ROOT=$star_root"; fi
  if [ "$dry_run" = "1" ]; then dry_run_arg="DRY_RUN=1"; fi
  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_splicing_merge_star_sj CONFIG="$CONFIG" SAMPLE="$sample" $star_root_arg $dry_run_arg $force_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_splicing_merge_star_sj CONFIG="$CONFIG" $star_root_arg $dry_run_arg $force_arg
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

run_cleanup() {
  local sample="${1:-}"
  local execute="${2:-0}"
  local threads="${3:-}"
  local nodes="${4:-}"
  local ppn="${5:-}"
  local mem="${6:-}"
  local walltime="${7:-}"
  local execute_arg=""
  local threads_arg=""
  local nodes_arg=""
  local ppn_arg=""
  local mem_arg=""
  local walltime_arg=""
  if [ "$execute" = "1" ]; then execute_arg="EXECUTE=1"; fi
  if [ -n "$threads" ]; then threads_arg="THREADS=$threads"; fi
  if [ -n "$nodes" ]; then nodes_arg="CLEANUP_NODES=$nodes"; fi
  if [ -n "$ppn" ]; then ppn_arg="CLEANUP_PPN=$ppn"; fi
  if [ -n "$mem" ]; then mem_arg="CLEANUP_MEM=$mem"; fi
  if [ -n "$walltime" ]; then walltime_arg="CLEANUP_WALLTIME=$walltime"; fi
  if [ -n "$sample" ]; then
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_cleanup_pre_mupexi CONFIG="$CONFIG" SAMPLE="$sample" $execute_arg $threads_arg $nodes_arg $ppn_arg $mem_arg $walltime_arg $force_arg $skip_running_arg
  else
    PIPELINE_DEFAULTS="$PIPELINE_DEFAULTS" make -C "$REPO" run_cleanup_pre_mupexi CONFIG="$CONFIG" $execute_arg $threads_arg $nodes_arg $ppn_arg $mem_arg $walltime_arg $force_arg $skip_running_arg
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
  local out_normal out_rna outdir germdir rna5_stranded_ext
  out_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  out_rna="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  rna5_stranded_ext="${rna5_stranded_vcf_extension:-rna5.1.filtered_edit_labeled.knownsites_summarised_qced_stranded.vcf.gz}"
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
    rna5.1) printf '%s\n' "${outdir}/${patient}_${rna5_stranded_ext}" ;;
    rna6) printf '%s\n' "${outdir}/${patient}_${rna6_merged_vcf_extension}" ;;
    rna7.0)
      rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
      printf '%s\n' "${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna7_smfixed_bam_suffix}"
      ;;
    rna7) printf '%s\n' "${outdir}/${patient}_${rna7_phased_vcf_extension}" ;;
    rna7.1) printf '%s\n' "${outdir}/${patient}_${rna7_phased_vcf_extension%.vcf.gz}.rna7.1.strand_filter_stats.tsv" ;;
    *) return 1 ;;
  esac
}

step_input_status() {
  local patient="$1"
  local step="$2"
  local out_normal out_rna dna_label rna_dir_label outdir germdir rna5_stranded_ext
  local rna_vcf_dir dna_vcf_dir srna_pref sdna_pref srna sdna
  out_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  out_rna="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
  rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
  rna5_stranded_ext="${rna5_stranded_vcf_extension:-rna5.1.filtered_edit_labeled.knownsites_summarised_qced_stranded.vcf.gz}"
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
    rna5.1)
      local in_rna51_vcf="${outdir}/${patient}_${rna5_qced_vcf_extension}"
      local in_rna51_bam="${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna7_smfixed_bam_suffix}"
      if [ -f "$in_rna51_vcf" ] && [ -s "$in_rna51_vcf" ] && [ -f "${GTF:-}" ] && [ -f "$in_rna51_bam" ] && [ -s "$in_rna51_bam" ]; then
        printf 'INPUT_OK\tVCF:%s ; BAM:%s ; GTF:%s\n' "$in_rna51_vcf" "$in_rna51_bam" "$GTF"
      elif [ ! -f "$in_rna51_vcf" ] || [ ! -s "$in_rna51_vcf" ]; then
        printf 'NO_INPUT\tVCF:%s\n' "$in_rna51_vcf"
      elif [ -z "${GTF:-}" ] || [ ! -f "$GTF" ]; then
        printf 'NO_INPUT\tGTF:%s\n' "${GTF:-CONFIG:GTF missing}"
      else
        printf 'NO_INPUT\tBAM:%s\n' "$in_rna51_bam"
      fi
      ;;
    rna6)
      local in_rna6="${outdir}/${patient}_${rna5_stranded_ext}"
      if [ ! -f "$in_rna6" ] || [ ! -s "$in_rna6" ]; then
        in_rna6="${outdir}/${patient}_${rna5_qced_vcf_extension}"
      fi
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
    rna7.1)
      local in_rna71_vcf="${outdir}/${patient}_${rna7_phased_vcf_extension}"
      local in_rna71_backup="${outdir}/${patient}_${rna7_phased_vcf_extension%.vcf.gz}.pre_rna7.1.vcf.gz"
      local in_rna71_bam="${bamdir}/${patient}_${rna_tumor_label:-RNA_TUMOR}/${patient}_${rna7_smfixed_bam_suffix}"
      local in_rna71_source="${in_rna71_backup}"
      if [ ! -f "$in_rna71_source" ] || [ ! -s "$in_rna71_source" ]; then
        in_rna71_source="${in_rna71_vcf}"
      fi
      if [ -f "$in_rna71_source" ] && [ -s "$in_rna71_source" ] && [ -f "$in_rna71_bam" ] && [ -s "$in_rna71_bam" ] && [ -f "${GTF:-}" ]; then
        printf 'INPUT_OK\tVCF:%s ; BAM:%s ; GTF:%s\n' "$in_rna71_source" "$in_rna71_bam" "$GTF"
      elif [ ! -f "$in_rna71_source" ] || [ ! -s "$in_rna71_source" ]; then
        printf 'NO_INPUT\tVCF:%s\n' "$in_rna71_source"
      elif [ -z "${GTF:-}" ] || [ ! -f "$GTF" ]; then
        printf 'NO_INPUT\tGTF:%s\n' "${GTF:-CONFIG:GTF missing}"
      else
        printf 'NO_INPUT\tBAM:%s\n' "$in_rna71_bam"
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
    rna5.1) printf '%s\n' "rna5.1_FilterByStrandness" ;;
    rna6) printf '%s\n' "rna6_MergeDnaRnaVcfs" ;;
    rna7.0) printf '%s\n' "rna7.0_FixRnaBamReadGroups" ;;
    rna7) printf '%s\n' "rna7_GenotypeAndPhaseMergedVcf" ;;
    rna7.1) printf '%s\n' "rna7.1_FilterByStrandness" ;;
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
    rna1|rna2|rna3|rna4|rna5|rna5.1|rna6|rna7.0|rna7|rna7.1|gdna1|gdna2|gdna3|gdna4) ;;
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
    rna1|rna2|rna3|rna4|rna5|rna5.1|rna6|rna7.0|rna7|rna7.1|gdna1|gdna2|gdna3|gdna4) ;;
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
  splicing)
    task="${1:-}"
    shift || true
    case "$task" in
      merge-star-sj)
        sample=""
        star_root=""
        dry_run="0"
        if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
          sample="$1"
          shift
        fi
        while [ $# -gt 0 ]; do
          case "${1:-}" in
            --root) star_root="${2:-}"; shift 2 ;;
            --dry-run) dry_run="1"; shift ;;
            *) echo "Unknown splicing option: $1" >&2; exit 1 ;;
          esac
        done
        run_splicing_merge_star_sj "$sample" "$star_root" "$dry_run"
        ;;
      *)
        echo "Unknown splicing task: ${task:-<missing>}" >&2
        exit 1
        ;;
    esac
    ;;
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
  cleanup)
    sample=""
    cleanup_execute="0"
    cleanup_threads=""
    cleanup_nodes=""
    cleanup_ppn=""
    cleanup_mem=""
    cleanup_walltime=""
    if [ $# -gt 0 ] && [[ "${1:-}" != -* ]]; then
      sample="$1"
      shift
    fi
    while [ $# -gt 0 ]; do
      case "${1:-}" in
        --execute) cleanup_execute="1"; shift ;;
        --threads) cleanup_threads="${2:-}"; shift 2 ;;
        --nodes) cleanup_nodes="${2:-}"; shift 2 ;;
        --ppn) cleanup_ppn="${2:-}"; shift 2 ;;
        --mem) cleanup_mem="${2:-}"; shift 2 ;;
        --walltime) cleanup_walltime="${2:-}"; shift 2 ;;
        *) echo "Unknown cleanup option: $1" >&2; exit 1 ;;
      esac
    done
    run_cleanup "$sample" "$cleanup_execute" "$cleanup_threads" "$cleanup_nodes" "$cleanup_ppn" "$cleanup_mem" "$cleanup_walltime"
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
      vep-dedup|vep-summary)
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
        run_research_vep_dedup "$sample" "$outdir"
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
      rna1|rna2|rna3|rna4|rna5|rna5.1|rna6|rna7.0|rna7|rna7.1|gdna1|gdna2|gdna3|gdna4)
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
