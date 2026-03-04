#!/usr/bin/bash
set -euo pipefail

# 4.1 Only RNA VCF: remove variants from RNA whose CHROM+POS exist in DNA tumour VCF.

usage() {
  cat <<'USAGE'
Usage: bash 4.1_OnlyRnaVcf.sh -c CONFIG [-s SAMPLE] [-f]
USAGE
}

force=0
while :; do
  case ${1:-} in
    -c|--config)
      if [ -n "${2:-}" ]; then config=$2; shift; else echo "ERROR: -c/--config requires a path" >&2; exit 1; fi
      ;;
    -s|--sample)
      if [ -n "${2:-}" ]; then sample=$2; shift; else echo "ERROR: -s/--sample requires a sample name" >&2; exit 1; fi
      ;;
    -f|--force)
      force=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      break
      ;;
  esac
  shift
done

[ -n "${config:-}" ] || { echo "ERROR: Config file needed. Use -c </path/to/CONFIG>." >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: Cannot find config file: $config" >&2; exit 1; }
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${rna_only_vcf_extension:?CONFIG must define rna_only_vcf_extension}"

sample_is_requested() {
  local sample_id="$1"
  local patient_id="$2"
  local requested="${sample:-}"
  [ -z "$requested" ] || [ "$sample_id" = "$requested" ] || [ "$patient_id" = "$requested" ]
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
    "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do
    value="${value%_${label}}"
  done
  printf '%s\n' "$value"
}

if [ -z "${sample:-}" ]; then
  echo "Running 4.1 for all samples in $samples"
else
  echo "Running 4.1 only for $sample"
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
  case "$line" in
    [[:space:]]*'#'*) continue ;;
  esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,	 ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1
  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"

  rna_vcf_dir="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"
  dna_vcf_dir="${vcfdir}/${name}_${dna_label}_vs_${name}_${out_normal_label}"

  sdna_vcf="${dna_vcf_dir}/${name}_${dna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz"
  srna_vcf="${rna_vcf_dir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz"

  rna_only_vcf="${rna_vcf_dir}/${name}_${rna_only_vcf_extension}"
  rna_only_log="${rna_only_vcf%.vcf.gz}.log.txt"

  if [ "$force" -eq 0 ] && [ -f "$rna_only_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $rna_only_vcf (use -f to overwrite)"
    continue
  fi
  job_name="${prefix}.${name}"
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
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load ngs tools htslib/1.23 bcftools/1.23 gatk/4.4.0.0
SCRIPT
    printf 'rna_vcf_dir=%q\n' "$rna_vcf_dir"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'srna_vcf=%q\n' "$srna_vcf"
    printf 'sdna_vcf=%q\n' "$sdna_vcf"
    printf 'rna_only_vcf=%q\n' "$rna_only_vcf"
    printf 'rna_only_log=%q\n' "$rna_only_log"
    cat <<'SCRIPT'
mkdir -p "$rna_vcf_dir"

bash "$rnae_scripts/rnae1_filter_somatic_from_rnavcf.sh" \
  --rna "$srna_vcf" --dna "$sdna_vcf" --out "$rna_only_vcf" 2>&1 | tee "$rna_only_log"

bcftools index -t "$rna_only_vcf"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=4,mem=8gb,walltime="00:08:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
