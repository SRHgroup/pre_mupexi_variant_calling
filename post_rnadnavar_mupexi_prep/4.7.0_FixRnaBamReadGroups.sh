#!/usr/bin/bash
set -euo pipefail

# 4.7.0 Optional: fix RNA BAM read group SM tags.

usage() {
  cat <<'USAGE'
Usage: bash 4.7.0_FixRnaBamReadGroups.sh -c CONFIG [-s SAMPLE] [-f]
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

: "${samples:?CONFIG must define samples}"
: "${bamdir:?CONFIG must define bamdir}"
: "${rna_bam_smfixed_suffix:?CONFIG must define rna_bam_smfixed_suffix}"

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

if [ -z "${sample:-}" ]; then
  echo "Running 4.7.0 for all samples in $samples"
else
  echo "Running 4.7.0 only for $sample"
fi

prefix=$(basename "${BASH_SOURCE[0]}" .sh)
scriptdir="${bamdir}/${prefix}.logs_and_reports"
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
  rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"

  inbam="${bamdir}/${name}_${rna_dir_label}/${name}_${rna_dir_label}.md.bam"
  outbam="${bamdir}/${name}_${rna_dir_label}/${name}_${rna_bam_smfixed_suffix}"
  header_sam="${outbam}.header.sam"
  check_txt="${outbam}.sm_check.txt"

  if [ "$force" -eq 0 ] && [ -f "$outbam" ]; then
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load tools htslib/1.23 samtools/1.23
SCRIPT
    printf 'inbam=%q\n' "$inbam"
    printf 'outbam=%q\n' "$outbam"
    printf 'header_sam=%q\n' "$header_sam"
    printf 'check_txt=%q\n' "$check_txt"
    printf 'rna_dir_label=%q\n' "$rna_dir_label"
    printf 'out_rna_label=%q\n' "$out_rna_label"
    cat <<'SCRIPT'
if [ ! -f "$inbam" ]; then
  echo "ERROR: missing input BAM: $inbam" >&2
  exit 1
fi

mkdir -p "$(dirname "$outbam")"

escaped_rna_label=$(printf '%s' "$rna_dir_label" | sed -e 's/[][(){}.^$*+?|/]/\\&/g')
samtools view -H "$inbam" | sed -E "s/SM:[^\t]*_?(${escaped_rna_label})(\\.[0-9]{4})?/SM:${out_rna_label}/g" > "$header_sam"
samtools reheader "$header_sam" "$inbam" > "$outbam"
samtools index -@ 8 "$outbam"

samtools view -H "$outbam" | awk -F'\t' '$1=="@RG"{for(i=1;i<=NF;i++) if($i~/^SM:/){sub(/^SM:/,"",$i); print $i}}' | sort -u | tee "$check_txt"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=8,mem=24gb,walltime="00:08:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
