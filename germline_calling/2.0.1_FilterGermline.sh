#!/usr/bin/bash
set -euo pipefail

# 2.0.1 Filter germline calls by QC thresholds.

usage() {
  cat <<'USAGE'
Usage: bash 2.0.1_FilterGermline.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
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
: "${output_extension_20:?CONFIG must define output_extension_20}"
: "${output_extension_201:?CONFIG must define output_extension_201}"

min_dp="${germline_min_dp:-10}"
min_qd="${germline_min_qd:-2.0}"
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
  germvcf="${germline_dir}/${name}_${output_extension_20}"
  filteredvcf="${germline_dir}/${name}_${output_extension_201}"

  if [ "$force" -eq 0 ] && [ -f "$filteredvcf" ]; then
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load tools ngs java/17-openjdk gatk/4.4.0.0
SCRIPT
    printf 'germvcf=%q\n' "$germvcf"
    printf 'filteredvcf=%q\n' "$filteredvcf"
    printf 'min_dp=%q\n' "$min_dp"
    printf 'min_qd=%q\n' "$min_qd"
    printf 'germline_dir=%q\n' "$germline_dir"
    cat <<'SCRIPT'
if [ ! -f "$germvcf" ]; then
  echo "ERROR: missing germline VCF from 2.0: $germvcf" >&2
  exit 1
fi

gatk VariantFiltration \
  -V "$germvcf" \
  -O "$filteredvcf" \
  --filter-expression "DP < ${min_dp}" --filter-name "LowDepth" \
  --filter-expression "QD < ${min_qd}" --filter-name "LowQualityByDepth" \
  --cluster-size 3 --cluster-window-size 15
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=8,mem=16gb,walltime="00:06:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
