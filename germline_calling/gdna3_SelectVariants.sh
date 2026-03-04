#!/usr/bin/bash
set -euo pipefail

# gdna3: keep only PASS/uncalled-filter variants from filtered germline VCF.

usage() {
  cat <<'USAGE'
Usage: bash gdna3_SelectVariants.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
USAGE
}

force=0
while :; do
  case ${1:-} in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c/--config requires a path" >&2; exit 1; }
      config=$2; shift ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s/--sample requires a value" >&2; exit 1; }
      sample=$2; shift ;;
    -f|--force) force=1 ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
  shift
done

[ -n "${config:-}" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${gdna2_vcf_extension:?CONFIG must define gdna2_vcf_extension}"
: "${gdna3_vcf_extension:?CONFIG must define gdna3_vcf_extension}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label:-DNA_NORMAL}" "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}"
    "${out_dna_normal_label:-DNA_NORMAL}" "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
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

  germline_dir="${vcfdir}/${name}_${out_normal_label}"
  filteredvcf="${germline_dir}/${name}_${gdna2_vcf_extension}"
  selectedvcf="${germline_dir}/${name}_${gdna3_vcf_extension}"

  if [ "$force" -eq 0 ] && [ -f "$selectedvcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $selectedvcf (use -f to overwrite)"
    continue
  fi
  if [ ! -s "$filteredvcf" ]; then
    echo "[precheck] ${prefix}.${name}: missing/empty filtered germline VCF from gdna2: $filteredvcf" >&2
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
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
module load tools ngs java/17-openjdk gatk/4.4.0.0
SCRIPT
    printf 'filteredvcf=%q\n' "$filteredvcf"
    printf 'selectedvcf=%q\n' "$selectedvcf"
    printf 'germline_dir=%q\n' "$germline_dir"
    cat <<'SCRIPT'
if [ ! -f "$filteredvcf" ]; then
  echo "ERROR: missing filtered germline VCF from gdna2: $filteredvcf" >&2
  exit 1
fi

gatk SelectVariants \
  -V "$filteredvcf" \
  -O "$selectedvcf" \
  --exclude-filtered true
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=8,mem=16gb,walltime="00:04:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
