#!/usr/bin/bash
set -euo pipefail

# rna4: collapse multiple RNA tumour samples into one RNA sample.

usage() {
  cat <<'USAGE'
Usage: bash rna4_SummariseRnaMetrics.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${rna3_knownsites_vcf_extension:?CONFIG must define rna3_knownsites_vcf_extension}"
: "${rna4_summarised_vcf_extension:?CONFIG must define rna4_summarised_vcf_extension}"

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

precheck_vcfgz() {
  local path="$1"
  local tag="$2"
  local hint="$3"
  if [ ! -f "$path" ]; then
    echo "[precheck] ${tag}: missing input VCF: $path (${hint})" >&2
    return 1
  fi
  if ! gzip -t "$path" 2>/dev/null; then
    echo "[precheck] ${tag}: invalid gzip VCF: $path" >&2
    return 1
  fi
  if zgrep -m1 '^#CHROM[[:space:]]' "$path" >/dev/null 2>&1; then
    return 0
  fi
  local first_non_meta
  first_non_meta="$(gzip -dc "$path" 2>/dev/null | awk '
    BEGIN{found=0}
    /^##/ {next}
    {
      sub(/\r$/,"",$0)
      sub(/^\xef\xbb\xbf/,"",$0)
      print
      found=1
      exit
    }
    END{ if (!found) exit 2 }
  ' || true)"
  if [[ "$first_non_meta" =~ ^#CHROM[[:space:]] ]]; then
    return 0
  fi
  if ! zgrep -m1 '^#CHROM' "$path" >/dev/null 2>&1; then
    echo "[precheck] ${tag}: malformed VCF header (missing #CHROM): $path" >&2
    return 1
  fi
}

# Use shared precheck implementation everywhere.
# shellcheck disable=SC1090
source "${repo_root}/post_rnadnavar_mupexi_prep/lib/vcf_precheck.sh"

if [ -z "${sample:-}" ]; then
  echo "Running rna4 for all samples in $samples"
else
  echo "Running rna4 only for $sample"
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
  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  annot_vcf="${outdir_only}/${name}_${rna3_knownsites_vcf_extension}"
  out_vcf="${outdir_only}/${name}_${rna4_summarised_vcf_extension}"
  stats_file="${out_vcf%.vcf.gz}.samples_and_counts.tsv"

  if [ "$force" -eq 0 ] && [ -f "$out_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $out_vcf (use -f to overwrite)"
    continue
  fi
  if ! precheck_vcfgz "$annot_vcf" "${prefix}.${name}" "run rna3 first"; then
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
    printf 'annot_vcf=%q\n' "$annot_vcf"
    printf 'out_vcf=%q\n' "$out_vcf"
    printf 'stats_file=%q\n' "$stats_file"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'dna_normal_label=%q\n' "${dna_normal_label:-DNA_NORMAL}"
    printf 'rna_tumor_label=%q\n' "${rna_tumor_label:-RNA_TUMOR}"
    printf 'out_dna_normal_label=%q\n' "${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
    printf 'out_rna_tumor_label=%q\n' "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    cat <<'SCRIPT'
if [ ! -f "$annot_vcf" ]; then
  echo "ERROR: missing input (run rna3 first): $annot_vcf" >&2
  exit 1
fi

python3 "$rnae_scripts/rnae4_summarise_RNA_metrics.py" \
  -i "$annot_vcf" -o "$out_vcf" \
  --normal-label "$dna_normal_label" \
  --rna-label "$rna_tumor_label" \
  --out-normal-name "$out_dna_normal_label" \
  --out-rna-name "$out_rna_tumor_label"

if ! bcftools view -h "$out_vcf" | grep -q '^#CHROM'; then
  echo "ERROR: malformed output VCF from 4.4 (missing #CHROM): $out_vcf" >&2
  rm -f "$out_vcf" "$out_vcf.tbi"
  exit 1
fi

bcftools index -t "$out_vcf"

in_total=$(bcftools view -H "$annot_vcf" | wc -l)
out_total=$(bcftools view -H "$out_vcf" | wc -l)

{
  printf "metric\tvalue\n"
  printf "input_total_variants\t%s\n" "$in_total"
  printf "output_total_variants\t%s\n" "$out_total"
  printf "input_samples\t%s\n" "$(bcftools query -l "$annot_vcf" | tr '\n' ' ')"
  printf "output_samples\t%s\n" "$(bcftools query -l "$out_vcf" | tr '\n' ' ')"
} > "$stats_file"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=4,mem=16gb,walltime="00:08:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
