#!/usr/bin/bash
set -euo pipefail

# 3.0 Filter selected germline calls by adjacency to somatic variants (optional RNA variants).

usage() {
  cat <<'USAGE'
Usage: bash 3.0_FilterGermlineByAdjacency.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
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
: "${output_extension_202:?CONFIG must define output_extension_202}"
: "${output_extension_30:?CONFIG must define output_extension_30}"

distance="${germline_proximity_distance:-33}"
include_rna="${germline_include_rna_proximity:-0}"

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

helpers_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"

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
  germ_in="${germline_dir}/${name}_${output_extension_202}"
  germ_out="${germline_dir}/${name}_${output_extension_30}"

  dna_dir="${vcfdir}/${name}_${out_dna_label}_vs_${name}_${out_normal_label}"
  somatic_vcf="${dna_dir}/${name}_${out_dna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz"

  rna_vcf=""
  if [ "$include_rna" = "1" ]; then
    rna_dir="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"
    rna_vcf="${rna_dir}/${name}_${rna_vcf_knownsites_extension}"
  fi

  if [ "$force" -eq 0 ] && [ -f "$germ_out" ]; then
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load tools ngs anaconda3/2025.06-1
SCRIPT
    printf 'germ_in=%q\n' "$germ_in"
    printf 'germ_out=%q\n' "$germ_out"
    printf 'germline_dir=%q\n' "$germline_dir"
    printf 'somatic_vcf=%q\n' "$somatic_vcf"
    printf 'rna_vcf=%q\n' "$rna_vcf"
    printf 'distance=%q\n' "$distance"
    printf 'helpers_dir=%q\n' "$helpers_dir"
    cat <<'SCRIPT'
if [ ! -f "$germ_in" ]; then
  echo "ERROR: missing germline VCF from 2.0.2: $germ_in" >&2
  exit 1
fi
if [ ! -f "$somatic_vcf" ]; then
  echo "ERROR: missing DNA somatic VCF for proximity filter: $somatic_vcf" >&2
  exit 1
fi

if [ -n "$rna_vcf" ] && [ -f "$rna_vcf" ]; then
  python3 "$helpers_dir/filter_germline_by_proximity.py" -g "$germ_in" -s "$somatic_vcf" --rna "$rna_vcf" -d "$distance" -o "$germ_out"
else
  python3 "$helpers_dir/filter_germline_by_proximity.py" -g "$germ_in" -s "$somatic_vcf" -d "$distance" -o "$germ_out"
fi
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=4,mem=8gb,walltime="00:04:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
