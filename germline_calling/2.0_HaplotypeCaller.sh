#!/usr/bin/bash
set -euo pipefail

# 2.0 Call germline variants in DNA normal BAM with GATK HaplotypeCaller.

usage() {
  cat <<'USAGE'
Usage: bash 2.0_HaplotypeCaller.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
USAGE
}

force=0
while :; do
  case ${1:-} in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c/--config requires a path" >&2; exit 1; }
      config=$2
      shift
      ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s/--sample requires a value" >&2; exit 1; }
      sample=$2
      shift
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

[ -n "${config:-}" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }
source "$config"

: "${samples:?CONFIG must define samples}"
: "${bamdir:?CONFIG must define bamdir}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${FASTA:?CONFIG must define FASTA}"
: "${output_extension_20:?CONFIG must define output_extension_20}"

dna_normal_label="${dna_normal_label:-DNA_NORMAL}"
dna_bam_suffix="${dna_bam_suffix:-md.bam}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label}}"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label}" "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}"
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

if [ -z "${sample:-}" ]; then
  echo "Running 2.0 for all patients in $samples"
else
  echo "Running 2.0 only for $sample"
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

  normal_bam="${bamdir}/${name}_${dna_normal_label}/${name}_${dna_normal_label}.${dna_bam_suffix}"
  germline_dir="${vcfdir}/${name}_${out_normal_label}"
  outvcf="${germline_dir}/${name}_${output_extension_20}"

  if [ "$force" -eq 0 ] && [ -f "$outvcf" ]; then
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load tools ngs java/17-openjdk gatk/4.4.0.0
SCRIPT
    printf 'normal_bam=%q\n' "$normal_bam"
    printf 'germline_dir=%q\n' "$germline_dir"
    printf 'outvcf=%q\n' "$outvcf"
    printf 'FASTA=%q\n' "$FASTA"
    cat <<'SCRIPT'
if [ ! -f "$normal_bam" ]; then
  echo "ERROR: missing DNA normal BAM: $normal_bam" >&2
  exit 1
fi

mkdir -p "$germline_dir"

gatk HaplotypeCaller \
  -I "$normal_bam" \
  -R "$FASTA" \
  -O "$outvcf" \
  -G StandardAnnotation
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=8,mem=24gb,walltime="00:12:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
