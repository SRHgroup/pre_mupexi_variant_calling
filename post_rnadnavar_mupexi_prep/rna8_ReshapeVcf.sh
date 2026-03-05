#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash rna8_ReshapeVcf.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${rna7_phased_vcf_extension:?CONFIG must define rna7_phased_vcf_extension}"

rna8_reshaped_vcf_extension="${rna8_reshaped_vcf_extension:-rna8.DNAt_DNAn_RNAt_merged_phased.reshaped.vcf.gz}"
rna8_tumor_output_label="${rna8_tumor_output_label:-TUMOR}"

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
  # NOTE:
  # Some valid cluster/bgzip VCFs have formatting quirks (BOM/whitespace/viewer issues)
  # that make simple grep-based #CHROM checks unreliable on login nodes.
  # Keep precheck strict only for existence + gzip integrity and let bcftools do the
  # authoritative parse inside the qsub job.
}

if [ -z "${sample:-}" ]; then
  echo "Running rna8 for all samples in $samples"
else
  echo "Running rna8 only for $sample"
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
  in_vcf="${outdir_only}/${name}_${rna7_phased_vcf_extension}"
  out_vcf="${outdir_only}/${name}_${rna8_reshaped_vcf_extension}"

  if [ "$force" -eq 0 ] && [ -f "$out_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $out_vcf (use -f to overwrite)"
    continue
  fi
  if ! precheck_vcfgz "$in_vcf" "${prefix}.${name}" "run rna7 first"; then
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
module load ${modules_rna:-ngs tools htslib/1.23 bcftools/1.23 anaconda3/2025.06-1}
SCRIPT
    printf 'in_vcf=%q\n' "$in_vcf"
    printf 'out_vcf=%q\n' "$out_vcf"
    printf 'out_normal=%q\n' "${out_dna_normal_label:-DNA_NORMAL}"
    printf 'out_dna=%q\n' "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
    printf 'out_rna=%q\n' "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    printf 'tumor_target=%q\n' "$rna8_tumor_output_label"
    cat <<'SCRIPT'
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# Find actual sample names in current VCF by suffix patterns.
normal_sample=$(bcftools query -l "$in_vcf" | awk -v n="$out_normal" 'toupper($0)==toupper(n){print; exit}')
dna_sample=$(bcftools query -l "$in_vcf" | awk -v d="$out_dna" 'toupper($0)==toupper(d){print; exit}')
rna_sample=$(bcftools query -l "$in_vcf" | awk -v r="$out_rna" 'toupper($0)==toupper(r){print; exit}')

[ -n "$normal_sample" ] || { echo "ERROR: normal sample not found in $in_vcf"; exit 1; }
[ -n "$dna_sample" ] || { echo "ERROR: DNA tumor sample not found in $in_vcf"; exit 1; }

bcftools view -s "${normal_sample},${dna_sample}" -Oz -o "${tmpdir}/subset.vcf.gz" "$in_vcf"
bcftools index -t "${tmpdir}/subset.vcf.gz"

# Rename selected samples to standardized output names.
{
  printf '%s\t%s\n' "$normal_sample" "$out_normal"
  printf '%s\t%s\n' "$dna_sample" "$tumor_target"
} > "${tmpdir}/rename.tsv"

bcftools reheader -s "${tmpdir}/rename.tsv" -o "${tmpdir}/reshaped.vcf.gz" "${tmpdir}/subset.vcf.gz"
bcftools index -t "${tmpdir}/reshaped.vcf.gz"
mv "${tmpdir}/reshaped.vcf.gz" "$out_vcf"
mv "${tmpdir}/reshaped.vcf.gz.tbi" "${out_vcf}.tbi"

echo "[done] wrote reshaped VCF: $out_vcf"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=4,mem=12gb,walltime="00:04:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"
  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
