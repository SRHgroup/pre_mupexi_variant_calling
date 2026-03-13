#!/usr/bin/bash
set -euo pipefail

# rna5.1: strand-aware filtering of RNA-edit VCF using transcript strand and RNA BAM.

usage() {
  cat <<'USAGE'
Usage: bash rna5.1_FilterByStrandness.sh -c CONFIG [-s SAMPLE] [-f]
USAGE
}

force=0
while :; do
  case ${1:-} in
    -c|--config)
      if [ -n "${2:-}" ]; then config=$2; shift; else echo "ERROR: -c/--config requires a path" >&2; exit 1; fi ;;
    -s|--sample)
      if [ -n "${2:-}" ]; then sample=$2; shift; else echo "ERROR: -s/--sample requires a sample name" >&2; exit 1; fi ;;
    -f|--force) force=1 ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
  shift
done

[ -n "${config:-}" ] || { echo "ERROR: Config file needed. Use -c </path/to/CONFIG>." >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: Cannot find config file: $config" >&2; exit 1; }
source "$config"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pipeline_defaults="${PIPELINE_DEFAULTS:-${repo_root}/pipeline_defaults/toolchain.defaults.sh}"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${bamdir:?CONFIG must define bamdir}"
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${GTF:?CONFIG must define GTF}"
: "${rna5_qced_vcf_extension:?CONFIG must define rna5_qced_vcf_extension}"

rna5_stranded_vcf_extension="${rna5_stranded_vcf_extension:-rna5.1.filtered_edit_labeled.knownsites_summarised_qced_stranded.vcf.gz}"
rna_strand_protocol="${rna_strand_protocol:-fr-firststrand}"
rna_strand_min_mapq="${rna_strand_min_mapq:-20}"
rna_strand_min_baseq="${rna_strand_min_baseq:-20}"
rna_strand_min_expected_frac="${rna_strand_min_expected_frac:-0.8}"

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
  local requested="${sample:-}"
  [ -z "$requested" ] || [ "$sample_id" = "$requested" ] || [ "$patient_id" = "$requested" ]
}

# Use shared precheck implementation everywhere.
# shellcheck disable=SC1090
source "${repo_root}/post_rnadnavar_mupexi_prep/lib/vcf_precheck.sh"

if [ -z "${sample:-}" ]; then
  echo "Running rna5.1 for all samples in $samples"
else
  echo "Running rna5.1 only for $sample"
fi

prefix=$(basename "${BASH_SOURCE[0]}" .sh)
scriptdir="${vcfdir}/${prefix}.logs_and_reports"
logdir="${scriptdir}/logs"
repdir="${scriptdir}/reports"
mkdir -p "$logdir" "$repdir"
qsub_depend_arg=()
if [ -n "${QSUB_DEPEND:-}" ]; then
  qsub_depend_arg=(-W "depend=${QSUB_DEPEND}")
fi

declare -A seen_patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,	 ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1

  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  in_vcf="${outdir_only}/${name}_${rna5_qced_vcf_extension}"
  out_vcf="${outdir_only}/${name}_${rna5_stranded_vcf_extension}"
  out_vcf_raw="${out_vcf%.gz}"
  support_tsv="${out_vcf%.vcf.gz}.strand_support.tsv"
  blacklist_tsv="${out_vcf%.vcf.gz}.strand_blacklist.tsv"
  stats_file="${out_vcf%.vcf.gz}.strand_filter_stats.tsv"
  rna_bam="${bamdir}/${name}_${rna_dir_label}/${name}_${rna7_smfixed_bam_suffix:-${rna_tumor_label:-RNA_TUMOR}.md.SMfixed.bam}"

  if [ "$force" -eq 0 ] && [ -f "$out_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $out_vcf (use -f to overwrite)"
    continue
  fi
  if ! precheck_vcfgz "$in_vcf" "${prefix}.${name}" "run rna5 first"; then
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi
  if [ ! -f "$rna_bam" ] || [ ! -s "$rna_bam" ]; then
    echo "[precheck] ${prefix}.${name}: missing RNA BAM: $rna_bam (run rna7.0 first)" >&2
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi
  if [ ! -f "$GTF" ]; then
    echo "[precheck] ${prefix}.${name}: missing GTF: $GTF" >&2
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi

  job_name="${prefix}.${name}"
  if [ "${SKIP_RUNNING:-0}" = "1" ]; then
    active_jobid=""
    if command -v qselect >/dev/null 2>&1; then
      active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "$job_name" 2>/dev/null | head -n1 || true)"
    fi
    if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
      active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="$job_name" '$4==n {print $1; exit}')"
    fi
    if [ -n "$active_jobid" ]; then
      echo "[skip] ${job_name}: job already active in scheduler: ${active_jobid}"
      continue
    fi
  fi
  submit_marker="${logdir}/submitted.${job_name}.jobid"
  if [ -f "$submit_marker" ]; then
    prev_jobid="$(head -n1 "$submit_marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ] && command -v qstat >/dev/null 2>&1 && qstat "$prev_jobid" >/dev/null 2>&1; then
      echo "[skip] ${job_name}: job already queued/running: ${prev_jobid}"
      continue
    fi
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    printf 'export PIPELINE_DEFAULTS=%q\n' "$pipeline_defaults"
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  source "$PIPELINE_DEFAULTS"
fi
module load ${modules_rna:-ngs tools htslib/1.23 bcftools/1.23 gatk/4.5.0.0 anaconda3/2025.06-1}
SCRIPT
    printf 'in_vcf=%q\n' "$in_vcf"
    printf 'out_vcf=%q\n' "$out_vcf"
    printf 'out_vcf_raw=%q\n' "$out_vcf_raw"
    printf 'support_tsv=%q\n' "$support_tsv"
    printf 'blacklist_tsv=%q\n' "$blacklist_tsv"
    printf 'stats_file=%q\n' "$stats_file"
    printf 'rna_bam=%q\n' "$rna_bam"
    printf 'gtf=%q\n' "$GTF"
    printf 'patient_name=%q\n' "$name"
    printf 'protocol=%q\n' "$rna_strand_protocol"
    printf 'min_mapq=%q\n' "$rna_strand_min_mapq"
    printf 'min_baseq=%q\n' "$rna_strand_min_baseq"
    printf 'min_expected_frac=%q\n' "$rna_strand_min_expected_frac"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    cat <<'SCRIPT'
if [ ! -f "$in_vcf" ]; then
  echo "ERROR: missing input (run rna5 first): $in_vcf" >&2
  exit 1
fi
if [ ! -f "$rna_bam" ]; then
  echo "ERROR: missing RNA BAM (run rna7.0 first): $rna_bam" >&2
  exit 1
fi
if [ ! -f "$gtf" ]; then
  echo "ERROR: missing GTF: $gtf" >&2
  exit 1
fi

python3 "$rnae_scripts/rnae5_1_filter_by_strandedness.py" \
  -i "$in_vcf" \
  -o "$out_vcf_raw" \
  --bam "$rna_bam" \
  --gtf "$gtf" \
  --patient "$patient_name" \
  --protocol "$protocol" \
  --min-mapq "$min_mapq" \
  --min-baseq "$min_baseq" \
  --min-expected-frac "$min_expected_frac" \
  --support-tsv "$support_tsv" \
  --blacklist-tsv "$blacklist_tsv" \
  --stats "$stats_file"

bgzip -f -c "$out_vcf_raw" > "$out_vcf"
rm -f "$out_vcf_raw"

if ! bcftools view -h "$out_vcf" | grep -q '^#CHROM'; then
  echo "ERROR: rna5.1 produced malformed VCF (missing #CHROM): $out_vcf" >&2
  rm -f "$out_vcf" "$out_vcf.tbi"
  exit 1
fi

bcftools index -t "$out_vcf"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=4,mem=20gb,walltime="00:12:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
