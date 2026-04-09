#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash bin/cleanup_pre_mupexi_outputs.sh -c CONFIG [-s PATIENT] [--execute] [--threads N]

Policy:
  - keep the first Mutect2 filtered RNA and DNA VCFs used by rna1
  - keep the final phased merged VCF used by MuPeXI (rna7_phased_vcf_extension)
  - if RNA final VCF is absent but dna_only_phased_vcf_extension exists, keep that final VCF instead
  - convert the last available pre-Mutect2 DNA normal, DNA tumor, and RNA tumor BAMs to CRAM
    (prefer recalibrated BAMs when present, otherwise fall back to markduplicates BAMs)
  - delete other heavy BAM/CRAM/VCF intermediates inside the per-patient BAM/VCF directories,
    including upstream recalibrated BAMs and STAR mapped BAM shards

Default mode is dry-run. Use --execute to make changes.

Options:
  -c, --config    Path to CONFIG file
  -s, --sample    Optional patient/sample base name to process
  --execute       Convert/delete files instead of printing the plan
  --threads       Threads for samtools CRAM conversion (default: 8)
  -h, --help      Show this help
USAGE
}

config=""
selected_sample=""
execute=0
threads=8

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c requires a path" >&2; exit 1; }
      config="$2"
      shift 2
      ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s requires a patient/sample" >&2; exit 1; }
      selected_sample="$2"
      shift 2
      ;;
    --execute)
      execute=1
      shift
      ;;
    --threads)
      [ -n "${2:-}" ] || { echo "ERROR: --threads requires a number" >&2; exit 1; }
      threads="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[ -n "$config" ] || { usage >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${bamdir:?CONFIG must define bamdir}"
[ -f "$samples" ] || { echo "ERROR: samples file not found: $samples" >&2; exit 1; }

if [ "$execute" -eq 1 ]; then
  : "${FASTA:?CONFIG must define FASTA for BAM->CRAM conversion}"
  if ! command -v samtools >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
      module load ${modules_rna_bamfix:-tools htslib/1.23 samtools/1.23} >/dev/null 2>&1 || true
    fi
  fi
  command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools is required for --execute" >&2; exit 1; }
fi

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
    "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do
    value="${value%_${label}}"
  done
  printf '%s\n' "$value"
}

sample_is_requested() {
  local sample_id="$1"
  local patient_id="$2"
  [ -z "$selected_sample" ] || [ "$selected_sample" = "$sample_id" ] || [ "$selected_sample" = "$patient_id" ]
}

resolve_patient_placeholder() {
  local template="$1"
  local patient="$2"
  printf '%s\n' "${template//\{patient\}/$patient}"
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
  local patient="$1"
  local sample_label="$2"
  local suffix="$3"
  local alt_label

  alt_label="$(flip_tumor_spelling "$sample_label")"
  printf '%s\n' \
    "${bamdir}/${patient}_${sample_label}/${patient}_${sample_label}.${suffix}" \
    "${bamdir}/${patient}_${sample_label}/${patient}_${suffix}" \
    "${bamdir}/${patient}_${sample_label}/${patient}_${sample_label}.md.bam" \
    "${bamdir}/${patient}_${sample_label}/${patient}_md.bam"

  if [ "$alt_label" != "$sample_label" ]; then
    printf '%s\n' \
      "${bamdir}/${patient}_${alt_label}/${patient}_${alt_label}.${suffix}" \
      "${bamdir}/${patient}_${alt_label}/${patient}_${suffix}" \
      "${bamdir}/${patient}_${alt_label}/${patient}_${alt_label}.md.bam" \
      "${bamdir}/${patient}_${alt_label}/${patient}_md.bam"
  fi
}

resolve_recal_bam_candidates() {
  local patient="$1"
  local sample_label="$2"
  local alt_label

  alt_label="$(flip_tumor_spelling "$sample_label")"
  printf '%s\n' \
    "${recaldir}/${patient}_${sample_label}/${patient}_${sample_label}.recal.bam" \
    "${recaldir}/${patient}_${sample_label}/${patient}.recal.bam"

  if [ "$alt_label" != "$sample_label" ]; then
    printf '%s\n' \
      "${recaldir}/${patient}_${alt_label}/${patient}_${alt_label}.recal.bam" \
      "${recaldir}/${patient}_${alt_label}/${patient}.recal.bam"
  fi
}

resolve_vcf_source_path() {
  local patient="$1"
  local kind="$2"
  local out_normal out_rna dna_label rna_vcf_dir dna_vcf_dir ext pref result

  out_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  out_rna="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
  rna_vcf_dir="${vcfdir}/${patient}_${out_rna}_vs_${patient}_${out_normal}"
  dna_vcf_dir="${vcfdir}/${patient}_${dna_label}_vs_${patient}_${out_normal}"

  if [ "$kind" = "rna" ]; then
    ext="${source_rna_mutect2_vcf_extension:-${out_rna}_vs_{patient}_${out_normal}.mutect2.filtered.vcf.gz}"
    ext="$(resolve_patient_placeholder "$ext" "$patient")"
    ext="${ext%\}}"
    pref="${rna_vcf_dir}/${patient}_${ext}"
    result="$(pick_first_existing \
      "$pref" \
      "${rna_vcf_dir}/${patient}_${out_rna}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      "${rna_vcf_dir}/${patient}_${rna_tumor_label:-RNA_TUMOR}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      "${rna_vcf_dir}/${patient}_RNA_TUMOR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      "${rna_vcf_dir}/${patient}_RNA_TUMOUR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
    )" || result=""
  else
    ext="${source_dna_mutect2_vcf_extension:-${dna_label}_vs_{patient}_${out_normal}.mutect2.filtered.vcf.gz}"
    ext="$(resolve_patient_placeholder "$ext" "$patient")"
    ext="${ext%\}}"
    pref="${dna_vcf_dir}/${patient}_${ext}"
    result="$(pick_first_existing \
      "$pref" \
      "${dna_vcf_dir}/${patient}_${dna_label}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      "${dna_vcf_dir}/${patient}_${dna_tumor_label:-DNA_TUMOR}_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      "${dna_vcf_dir}/${patient}_DNA_TUMOR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
      "${dna_vcf_dir}/${patient}_DNA_TUMOUR_vs_${patient}_${out_normal}.mutect2.filtered.vcf.gz" \
    )" || result=""
  fi

  printf '%s\n' "$result"
}

vcf_sidecars() {
  local vcf="$1"
  [ -n "$vcf" ] || return 0

  if [[ "$vcf" == *.vcf.gz ]]; then
    printf '%s\n' "${vcf}.tbi" "${vcf}.csi"
  elif [[ "$vcf" == *.vcf ]]; then
    printf '%s\n' "${vcf}.idx"
  fi
}

bam_index_candidates() {
  local bam="$1"
  printf '%s\n' "${bam}.bai" "${bam%.bam}.bai"
}

future_cram_path() {
  local bam="$1"
  if [[ "$bam" == *.bam ]]; then
    printf '%s\n' "${bam%.bam}.cram"
  else
    printf '%s\n' "${bam}.cram"
  fi
}

file_is_heavy_candidate() {
  local path="$1"
  case "$path" in
    *.bam|*.bai|*.cram|*.crai|*.vcf|*.vcf.gz|*.tbi|*.csi|*.idx|*.pileups.table|*.bamout|*.bamout.bai|*.f1r2.tar.gz|*.read-orientation-model.tar.gz)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

star_mapped_dirs() {
  local patient="$1"
  local sample_label="$2"
  local alt_label

  alt_label="$(flip_tumor_spelling "$sample_label")"
  printf '%s\n' \
    "${stardir}/${patient}/${patient}_${sample_label}/mapped" \
    "${stardir}/${patient}_${sample_label}/mapped"
  if [ "$alt_label" != "$sample_label" ]; then
    printf '%s\n' \
      "${stardir}/${patient}/${patient}_${alt_label}/mapped" \
      "${stardir}/${patient}_${alt_label}/mapped"
  fi
}

declare -a delete_paths=()
declare -a keep_paths=()

array_contains() {
  local needle="$1"
  shift || true
  local item
  for item in "$@"; do
    [ "$item" = "$needle" ] && return 0
  done
  return 1
}

reset_patient_state() {
  keep_paths=()
  delete_paths=()
}

add_keep_path() {
  local path="$1"
  [ -n "$path" ] || return 0
  if ! array_contains "$path" "${keep_paths[@]:-}"; then
    keep_paths+=("$path")
  fi
}

add_delete_path() {
  local path="$1"
  [ -n "$path" ] || return 0
  if ! array_contains "$path" "${delete_paths[@]:-}"; then
    delete_paths+=("$path")
  fi
}

gather_delete_candidates_from_dirs() {
  local dir path
  local -a dirs=("$@")

  for dir in "${dirs[@]}"; do
    [ -d "$dir" ] || continue
    while IFS= read -r path; do
      file_is_heavy_candidate "$path" || continue
      array_contains "$path" "${keep_paths[@]:-}" && continue
      add_delete_path "$path"
    done < <(find "$dir" -maxdepth 1 -type f | sort)
  done
}

gather_star_delete_candidates() {
  local patient="$1"
  local sample_label="$2"
  local mapped_dir path

  while IFS= read -r mapped_dir; do
    [ -d "$mapped_dir" ] || continue
    while IFS= read -r path; do
      file_is_heavy_candidate "$path" || continue
      array_contains "$path" "${keep_paths[@]:-}" && continue
      add_delete_path "$path"
    done < <(find "$mapped_dir" -maxdepth 1 -type f | sort)
  done < <(star_mapped_dirs "$patient" "$sample_label")
}

convert_bam_to_cram() {
  local bam="$1"
  local cram tmp tmp_crai idx

  cram="$(future_cram_path "$bam")"
  tmp="${cram}.tmp.$$"
  tmp_crai="${tmp}.crai"

  if [ -f "$cram" ] && [ -s "$cram" ]; then
    samtools quickcheck "$cram"
  else
    echo "CONVERT $bam -> $cram"
    samtools view -@ "$threads" -T "$FASTA" -C -o "$tmp" "$bam"
    samtools index -@ "$threads" -o "$tmp_crai" "$tmp"
    samtools quickcheck "$tmp"
    mv "$tmp" "$cram"
    mv "$tmp_crai" "${cram}.crai"
  fi

  rm -f "$bam"
  while IFS= read -r idx; do
    rm -f "$idx"
  done < <(bam_index_candidates "$bam")
}

out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
dna_suffix="${dna_bam_suffix:-md.bam}"
preprocessing_root="$(dirname "$bamdir")"
recaldir="${recaldir:-${preprocessing_root}/recalibrated}"
stardir="${stardir:-${preprocessing_root}/star}"

patients_seen=0
patients_ready=0
patients_skipped=0
planned_conversions=0
planned_deletions=0
executed_conversions=0
executed_deletions=0
selected_matched=0

declare -a seen_patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
  patient="$(sample_base_name "$sample_name")"
  [ -n "$patient" ] || continue
  sample_is_requested "$sample_name" "$patient" || continue
  array_contains "$patient" "${seen_patients[@]:-}" && continue
  seen_patients+=("$patient")
  if [ -n "$selected_sample" ]; then
    selected_matched=1
  fi
  patients_seen=$((patients_seen + 1))

  reset_patient_state

  germ_dir="${vcfdir}/${patient}_${out_normal_label}"
  dna_vcf_dir="${vcfdir}/${patient}_${out_dna_label}_vs_${patient}_${out_normal_label}"
  rna_vcf_dir="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}"

  source_rna_vcf="$(resolve_vcf_source_path "$patient" "rna")"
  source_dna_vcf="$(resolve_vcf_source_path "$patient" "dna")"
  final_rna_vcf="${rna_vcf_dir}/${patient}_${rna7_phased_vcf_extension:-}"
  final_dna_only_vcf=""
  if [ -n "${dna_only_phased_vcf_extension:-}" ]; then
    final_dna_only_vcf="${dna_vcf_dir}/${patient}_${dna_only_phased_vcf_extension}"
  fi

  final_keep_vcf=""
  if [ -n "${rna7_phased_vcf_extension:-}" ] && [ -f "$final_rna_vcf" ]; then
    final_keep_vcf="$final_rna_vcf"
  elif [ -n "$final_dna_only_vcf" ] && [ -f "$final_dna_only_vcf" ]; then
    final_keep_vcf="$final_dna_only_vcf"
  fi

  dna_normal_md_bam="$(pick_first_existing $(resolve_bam_candidates "$patient" "$out_normal_label" "$dna_suffix"))" || dna_normal_md_bam=""
  dna_tumor_md_bam="$(pick_first_existing $(resolve_bam_candidates "$patient" "$out_dna_label" "$dna_suffix"))" || dna_tumor_md_bam=""
  rna_tumor_md_bam="$(pick_first_existing $(resolve_bam_candidates "$patient" "$rna_dir_label" "$dna_suffix"))" || rna_tumor_md_bam=""
  dna_normal_recal_bam="$(pick_first_existing $(resolve_recal_bam_candidates "$patient" "$out_normal_label"))" || dna_normal_recal_bam=""
  dna_tumor_recal_bam="$(pick_first_existing $(resolve_recal_bam_candidates "$patient" "$out_dna_label"))" || dna_tumor_recal_bam=""
  rna_tumor_recal_bam="$(pick_first_existing $(resolve_recal_bam_candidates "$patient" "$rna_dir_label"))" || rna_tumor_recal_bam=""

  dna_normal_bam="${dna_normal_recal_bam:-$dna_normal_md_bam}"
  dna_tumor_bam="${dna_tumor_recal_bam:-$dna_tumor_md_bam}"
  rna_tumor_bam="${rna_tumor_recal_bam:-$rna_tumor_md_bam}"

  missing_reason=""
  if [ -z "$source_rna_vcf" ]; then
    missing_reason="missing source RNA Mutect2 filtered VCF"
  elif [ -z "$source_dna_vcf" ]; then
    missing_reason="missing source DNA Mutect2 filtered VCF"
  elif [ -z "$final_keep_vcf" ]; then
    missing_reason="missing final phased merged VCF"
  elif [ -z "$dna_normal_bam" ]; then
    missing_reason="missing DNA normal BAM"
  elif [ -z "$dna_tumor_bam" ]; then
    missing_reason="missing DNA tumor BAM"
  elif [ -z "$rna_tumor_bam" ]; then
    missing_reason="missing RNA tumor BAM"
  fi

  if [ -n "$missing_reason" ]; then
    patients_skipped=$((patients_skipped + 1))
    echo "SKIP $patient: $missing_reason"
    continue
  fi

  add_keep_path "$source_rna_vcf"
  add_keep_path "$source_dna_vcf"
  add_keep_path "$final_keep_vcf"

  while IFS= read -r sidecar; do
    [ -f "$sidecar" ] && add_keep_path "$sidecar"
  done < <(vcf_sidecars "$source_rna_vcf")
  while IFS= read -r sidecar; do
    [ -f "$sidecar" ] && add_keep_path "$sidecar"
  done < <(vcf_sidecars "$source_dna_vcf")
  while IFS= read -r sidecar; do
    [ -f "$sidecar" ] && add_keep_path "$sidecar"
  done < <(vcf_sidecars "$final_keep_vcf")

  dna_normal_cram="$(future_cram_path "$dna_normal_bam")"
  dna_tumor_cram="$(future_cram_path "$dna_tumor_bam")"
  rna_tumor_cram="$(future_cram_path "$rna_tumor_bam")"
  add_keep_path "$dna_normal_cram"
  add_keep_path "${dna_normal_cram}.crai"
  add_keep_path "$dna_tumor_cram"
  add_keep_path "${dna_tumor_cram}.crai"
  add_keep_path "$rna_tumor_cram"
  add_keep_path "${rna_tumor_cram}.crai"

  patient_dirs=()
  for dir in \
    "$germ_dir" \
    "$dna_vcf_dir" \
    "$rna_vcf_dir" \
    "$(dirname "${dna_normal_md_bam:-$dna_normal_bam}")" \
    "$(dirname "${dna_tumor_md_bam:-$dna_tumor_bam}")" \
    "$(dirname "${rna_tumor_md_bam:-$rna_tumor_bam}")" \
    "$(dirname "${dna_normal_recal_bam:-$dna_normal_bam}")" \
    "$(dirname "${dna_tumor_recal_bam:-$dna_tumor_bam}")" \
    "$(dirname "${rna_tumor_recal_bam:-$rna_tumor_bam}")" \
    "$(dirname "$dna_normal_bam")" \
    "$(dirname "$dna_tumor_bam")" \
    "$(dirname "$rna_tumor_bam")"
  do
    [ -d "$dir" ] || continue
    seen_dir=0
    for existing_dir in "${patient_dirs[@]:-}"; do
      if [ "$existing_dir" = "$dir" ]; then
        seen_dir=1
        break
      fi
    done
    [ "$seen_dir" -eq 1 ] || patient_dirs+=("$dir")
  done

  delete_paths=()
  gather_delete_candidates_from_dirs "${patient_dirs[@]}"
  gather_star_delete_candidates "$patient" "$rna_dir_label"

  patients_ready=$((patients_ready + 1))
  planned_conversions=$((planned_conversions + 3))
  planned_deletions=$((planned_deletions + ${#delete_paths[@]}))

  echo "PATIENT $patient"
  echo "  keep_vcfs:"
  echo "    $source_rna_vcf"
  echo "    $source_dna_vcf"
  echo "    $final_keep_vcf"
  echo "  keep_crams:"
  echo "    $dna_normal_cram"
  echo "    $dna_tumor_cram"
  echo "    $rna_tumor_cram"
  echo "  source_bams_for_cram:"
  echo "    $dna_normal_bam"
  echo "    $dna_tumor_bam"
  echo "    $rna_tumor_bam"
  echo "  delete_count: ${#delete_paths[@]}"

  if [ "$execute" -eq 0 ]; then
    for path in "${delete_paths[@]}"; do
      echo "    DELETE $path"
    done
    continue
  fi

  convert_bam_to_cram "$dna_normal_bam"
  executed_conversions=$((executed_conversions + 1))
  convert_bam_to_cram "$dna_tumor_bam"
  executed_conversions=$((executed_conversions + 1))
  convert_bam_to_cram "$rna_tumor_bam"
  executed_conversions=$((executed_conversions + 1))

  delete_paths=()
  gather_delete_candidates_from_dirs "${patient_dirs[@]}"
  gather_star_delete_candidates "$patient" "$rna_dir_label"
  for path in "${delete_paths[@]}"; do
    rm -f "$path"
    executed_deletions=$((executed_deletions + 1))
  done
done < "$samples"

echo
echo "SUMMARY"
echo "  mode: $([ "$execute" -eq 1 ] && echo execute || echo dry-run)"
echo "  patients_seen: $patients_seen"
echo "  patients_ready: $patients_ready"
echo "  patients_skipped: $patients_skipped"
echo "  planned_cram_conversions: $planned_conversions"
echo "  planned_deletions: $planned_deletions"
if [ "$execute" -eq 1 ]; then
  echo "  executed_cram_conversions: $executed_conversions"
  echo "  executed_deletions: $executed_deletions"
fi

if [ -n "$selected_sample" ]; then
  if [ "$selected_matched" -eq 0 ]; then
    echo "ERROR: selected patient/sample not found in SAMPLES: $selected_sample" >&2
    exit 3
  fi
  if [ "$patients_ready" -eq 0 ]; then
    exit 2
  fi
fi
