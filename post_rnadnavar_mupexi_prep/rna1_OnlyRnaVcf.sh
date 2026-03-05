#!/usr/bin/bash
set -euo pipefail

# rna1: remove variants from RNA whose CHROM+POS exist in DNA tumour VCF.

usage() {
  cat <<'USAGE'
Usage: bash rna1_OnlyRnaVcf.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${rna1_vcf_extension:?CONFIG must define rna1_vcf_extension}"

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

if [ -z "${sample:-}" ]; then
  echo "Running rna1 for all samples in $samples"
else
  echo "Running rna1 only for $sample"
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
  source_rna_mutect2_vcf_extension="${source_rna_mutect2_vcf_extension:-${out_rna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
  source_dna_mutect2_vcf_extension="${source_dna_mutect2_vcf_extension:-${dna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
  source_rna_mutect2_vcf_extension="$(resolve_patient_placeholder "$source_rna_mutect2_vcf_extension" "$name")"
  source_dna_mutect2_vcf_extension="$(resolve_patient_placeholder "$source_dna_mutect2_vcf_extension" "$name")"

  # Strip any accidental trailing braces from misconfigured templates.
  source_rna_mutect2_vcf_extension="${source_rna_mutect2_vcf_extension%\}}"
  source_dna_mutect2_vcf_extension="${source_dna_mutect2_vcf_extension%\}}"

  rna_vcf_dir="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"
  dna_vcf_dir="${vcfdir}/${name}_${dna_label}_vs_${name}_${out_normal_label}"

  srna_vcf_pref="${rna_vcf_dir}/${name}_${source_rna_mutect2_vcf_extension}"
  sdna_vcf_pref="${dna_vcf_dir}/${name}_${source_dna_mutect2_vcf_extension}"

  srna_vcf="$(pick_first_existing \
    "$srna_vcf_pref" \
    "${rna_vcf_dir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${name}_${rna_tumor_label:-RNA_TUMOR}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${name}_RNA_TUMOR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${name}_RNA_TUMOUR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
  )" || srna_vcf="$srna_vcf_pref"

  sdna_vcf="$(pick_first_existing \
    "$sdna_vcf_pref" \
    "${dna_vcf_dir}/${name}_${dna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${name}_${dna_tumor_label:-DNA_TUMOR}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${name}_DNA_TUMOR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${name}_DNA_TUMOUR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
  )" || sdna_vcf="$sdna_vcf_pref"

  rna_only_vcf="${rna_vcf_dir}/${name}_${rna1_vcf_extension}"
  rna_only_log="${rna_only_vcf%.vcf.gz}.log.txt"

  if [ "$force" -eq 0 ] && [ -f "$rna_only_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $rna_only_vcf (use -f to overwrite)"
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
