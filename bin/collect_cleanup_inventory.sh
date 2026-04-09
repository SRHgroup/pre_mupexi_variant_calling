#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash bin/collect_cleanup_inventory.sh -c CONFIG [-p PATIENT] [-n TOP_N]

Print a compact inventory for planning cleanup and BAM->CRAM conversion.

Options:
  -c, --config   Path to CONFIG file
  -p, --patient  Optional patient/sample base name to inspect
  -n, --top      Number of largest paths to show (default: 60)
  -h, --help     Show this help
USAGE
}

config=""
patient=""
top_n=60

while :; do
  case "${1:-}" in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c requires a path" >&2; exit 1; }
      config="$2"
      shift
      ;;
    -p|--patient)
      [ -n "${2:-}" ] || { echo "ERROR: -p requires a patient/sample name" >&2; exit 1; }
      patient="$2"
      shift
      ;;
    -n|--top)
      [ -n "${2:-}" ] || { echo "ERROR: -n requires a number" >&2; exit 1; }
      top_n="$2"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    "")
      break
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

[ -n "$config" ] || { usage >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${bamdir:?CONFIG must define bamdir}"
[ -f "$samples" ] || { echo "ERROR: samples file not found: $samples" >&2; exit 1; }

flip_tumor_spelling() {
  local value="$1"
  if [[ "$value" == *TUMOR* ]]; then
    printf '%s\n' "${value/TUMOR/TUMOUR}"
  elif [[ "$value" == *TUMOUR* ]]; then
    printf '%s\n' "${value/TUMOUR/TUMOR}"
  else
    printf '%s\n' "$value"
  fi
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

pick_first_patient() {
  local line sample_name base_name
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    case "$line" in
      [[:space:]]*'#'*) continue ;;
    esac
    sample_name=$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')
    base_name=$(sample_base_name "$sample_name")
    [ -n "$base_name" ] || continue
    printf '%s\n' "$base_name"
    return 0
  done < "$samples"
  return 1
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

resolve_bam_candidates() {
  local patient_name="$1"
  local sample_label="$2"
  local suffix="$3"
  local alt_label

  alt_label="$(flip_tumor_spelling "$sample_label")"
  printf '%s\n' \
    "${bamdir}/${patient_name}_${sample_label}/${patient_name}_${sample_label}.${suffix}" \
    "${bamdir}/${patient_name}_${sample_label}/${patient_name}_${suffix}" \
    "${bamdir}/${patient_name}_${sample_label}/${patient_name}_${sample_label}.md.bam" \
    "${bamdir}/${patient_name}_${sample_label}/${patient_name}_md.bam"

  if [ "$alt_label" != "$sample_label" ]; then
    printf '%s\n' \
      "${bamdir}/${patient_name}_${alt_label}/${patient_name}_${alt_label}.${suffix}" \
      "${bamdir}/${patient_name}_${alt_label}/${patient_name}_${suffix}" \
      "${bamdir}/${patient_name}_${alt_label}/${patient_name}_${alt_label}.md.bam" \
      "${bamdir}/${patient_name}_${alt_label}/${patient_name}_md.bam"
  fi
}

resolve_smfixed_candidates() {
  local patient_name="$1"
  local sample_label="$2"
  local suffix="${rna7_smfixed_bam_suffix:-${sample_label}.md.SMfixed.bam}"
  local alt_label alt_suffix

  alt_label="$(flip_tumor_spelling "$sample_label")"
  alt_suffix="$(flip_tumor_spelling "$suffix")"
  printf '%s\n' \
    "${bamdir}/${patient_name}_${sample_label}/${patient_name}_${suffix}" \
    "${bamdir}/${patient_name}_${sample_label}/${patient_name}_${alt_suffix}"

  if [ "$alt_label" != "$sample_label" ]; then
    printf '%s\n' \
      "${bamdir}/${patient_name}_${alt_label}/${patient_name}_${suffix}" \
      "${bamdir}/${patient_name}_${alt_label}/${patient_name}_${alt_suffix}"
  fi
}

print_section() {
  printf '\n== %s ==\n' "$1"
}

print_var() {
  local name="$1"
  local value="${!name-}"
  if [ -n "$value" ]; then
    printf '%s\t%s\n' "$name" "$value"
  else
    printf '%s\t%s\n' "$name" "(unset)"
  fi
}

print_path_status() {
  local label="$1"
  local path="$2"
  local size

  if [ -z "$path" ]; then
    printf 'unset\t%s\t-\t%s\n' "$label" "(unset)"
    return
  fi
  if [ -d "$path" ] || [ -f "$path" ]; then
    size=$(du -sh "$path" 2>/dev/null | awk '{print $1}')
    printf 'present\t%s\t%s\t%s\n' "$label" "${size:--}" "$path"
  else
    printf 'missing\t%s\t-\t%s\n' "$label" "$path"
  fi
}

patient="${patient:-$(pick_first_patient || true)}"
[ -n "$patient" ] || { echo "ERROR: could not determine a patient from $samples" >&2; exit 1; }

sample_lines=$(awk 'NF && $1 !~ /^#/ {count++} END {print count+0}' "$samples")

out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
dna_suffix="${dna_bam_suffix:-md.bam}"

germ_dir="${vcfdir}/${patient}_${out_normal_label}"
dna_vcf_dir="${vcfdir}/${patient}_${dna_label}_vs_${patient}_${out_normal_label}"
rna_vcf_dir="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}"

normal_bam="$(pick_first_existing $(resolve_bam_candidates "$patient" "$out_normal_label" "$dna_suffix"))" || normal_bam=""
tumor_bam="$(pick_first_existing $(resolve_bam_candidates "$patient" "$dna_label" "$dna_suffix"))" || tumor_bam=""
rna_md_bam="$(pick_first_existing $(resolve_bam_candidates "$patient" "$rna_dir_label" "$dna_suffix"))" || rna_md_bam=""
rna_smfixed_bam="$(pick_first_existing $(resolve_smfixed_candidates "$patient" "$rna_dir_label"))" || rna_smfixed_bam=""

normal_bam="${normal_bam:-${bamdir}/${patient}_${out_normal_label}/${patient}_${out_normal_label}.${dna_suffix}}"
tumor_bam="${tumor_bam:-${bamdir}/${patient}_${dna_label}/${patient}_${dna_label}.${dna_suffix}}"
rna_md_bam="${rna_md_bam:-${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna_dir_label}.${dna_suffix}}"
rna_smfixed_bam="${rna_smfixed_bam:-${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna7_smfixed_bam_suffix:-${rna_dir_label}.md.SMfixed.bam}}"

germ_vcf="${germ_dir}/${patient}_${gdna4_vcf_extension:-}"
rna_phased_vcf="${rna_vcf_dir}/${patient}_${rna7_phased_vcf_extension:-}"
dna_only_phased_vcf=""
if [ -n "${dna_only_phased_vcf_extension:-}" ]; then
  dna_only_phased_vcf="${dna_vcf_dir}/${patient}_${dna_only_phased_vcf_extension}"
fi

print_section "CONFIG SUMMARY"
print_var config
print_var samples
printf '%s\t%s\n' "sample_rows" "$sample_lines"
printf '%s\t%s\n' "selected_patient" "$patient"

print_section "KEY VARIABLES"
for name in \
  vcfdir bamdir FASTA dna_bam_suffix rna7_smfixed_bam_suffix gdna4_vcf_extension \
  rna7_phased_vcf_extension dna_only_phased_vcf_extension mupexi_outdir kaldir \
  mosdepthdir dna_normal_label dna_tumor_label rna_tumor_label out_dna_normal_label \
  out_dna_tumor_label out_rna_tumor_label
do
  print_var "$name"
done

print_section "ROOT SIZES"
print_path_status "bamdir" "$bamdir"
print_path_status "vcfdir" "$vcfdir"
print_path_status "kaldir" "${kaldir:-}"
print_path_status "mupexi_outdir" "${mupexi_outdir:-}"
print_path_status "mosdepthdir" "${mosdepthdir:-}"

roots=()
for root in "$bamdir" "$vcfdir" "${kaldir:-}" "${mupexi_outdir:-}" "${mosdepthdir:-}"; do
  [ -n "$root" ] || continue
  [ -e "$root" ] || continue
  roots+=("$root")
done

print_section "TOP SPACE HOGS"
if [ "${#roots[@]}" -eq 0 ]; then
  echo "No existing storage roots found from CONFIG."
else
  du -ah "${roots[@]}" 2>/dev/null | sort -hr | awk -v limit="$top_n" 'NR <= limit { print }'
fi

print_section "EXPECTED FINAL KEEPERS"
print_path_status "dna_normal_bam" "$normal_bam"
print_path_status "dna_tumor_bam" "$tumor_bam"
print_path_status "rna_md_bam" "$rna_md_bam"
print_path_status "rna7_smfixed_bam" "$rna_smfixed_bam"
print_path_status "gdna4_vcf" "$germ_vcf"
print_path_status "rna7_phased_vcf" "$rna_phased_vcf"
print_path_status "dna_only_phased_vcf" "$dna_only_phased_vcf"

print_section "PATIENT DIRECTORIES"
print_path_status "germline_dir" "$germ_dir"
print_path_status "dna_vcf_dir" "$dna_vcf_dir"
print_path_status "rna_vcf_dir" "$rna_vcf_dir"
print_path_status "dna_normal_dir" "$(dirname "$normal_bam")"
print_path_status "dna_tumor_dir" "$(dirname "$tumor_bam")"
print_path_status "rna_dir" "$(dirname "$rna_smfixed_bam")"

patient_dirs=()
add_patient_dir() {
  local dir="$1"
  local existing
  if [ ! -d "$dir" ]; then
    return 0
  fi
  for existing in "${patient_dirs[@]:-}"; do
    if [ "$existing" = "$dir" ]; then
      return 0
    fi
  done
  patient_dirs+=("$dir")
}

add_patient_dir "$germ_dir"
add_patient_dir "$dna_vcf_dir"
add_patient_dir "$rna_vcf_dir"
add_patient_dir "$(dirname "$normal_bam")"
add_patient_dir "$(dirname "$tumor_bam")"
add_patient_dir "$(dirname "$rna_md_bam")"
alt_normal_label="$(flip_tumor_spelling "$out_normal_label")"
alt_dna_label="$(flip_tumor_spelling "$dna_label")"
alt_rna_label="$(flip_tumor_spelling "$rna_dir_label")"
add_patient_dir "${bamdir}/${patient}_${alt_normal_label}"
add_patient_dir "${bamdir}/${patient}_${alt_dna_label}"
add_patient_dir "${bamdir}/${patient}_${alt_rna_label}"

print_section "PATIENT LAYOUT"
if [ "${#patient_dirs[@]}" -eq 0 ]; then
  echo "No patient-specific directories found for ${patient}."
else
  find "${patient_dirs[@]}" -maxdepth 2 \( -type d -o -type f \) | sort
fi
