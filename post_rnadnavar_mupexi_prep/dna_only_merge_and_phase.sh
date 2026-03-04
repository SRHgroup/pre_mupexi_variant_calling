#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash dna_only_merge_and_phase.sh -c CONFIG [-s SAMPLE] [-f] [--allow-with-rna]
USAGE
}

force=0
allow_with_rna=0
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
    --allow-with-rna)
      allow_with_rna=1
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
# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${bamdir:?CONFIG must define bamdir}"
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${FASTA:?CONFIG must define FASTA}"
: "${DICT:?CONFIG must define DICT}"

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
    "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
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

precheck_vcf_or_vcfgz() {
  local file="$1"
  local tag="$2"
  if [ ! -f "$file" ]; then
    echo "ERROR: missing input: $file" >&2
    echo "[skip] ${tag}: missing input" >&2
    return 1
  fi
  if [ ! -s "$file" ]; then
    echo "ERROR: empty input: $file" >&2
    echo "[skip] ${tag}: empty input" >&2
    return 1
  fi
  if [[ "$file" = *.gz ]]; then
    if ! gzip -t "$file" >/dev/null 2>&1; then
      echo "ERROR: corrupted/non-gzip input: $file" >&2
      echo "[skip] ${tag}: invalid input" >&2
      return 1
    fi
  fi
  return 0
}

if [ -z "${sample:-}" ]; then
  echo "Running dna-only merge+phase for all samples in $samples"
else
  echo "Running dna-only merge+phase only for $sample"
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

  sample_name=$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1

  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"

  rna_vcf_dir="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"
  dna_vcf_dir="${vcfdir}/${name}_${dna_label}_vs_${name}_${out_normal_label}"
  germ_dir="${vcfdir}/${name}_${out_normal_label}"

  source_rna_mutect2_vcf_extension="${source_rna_mutect2_vcf_extension:-${out_rna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
  source_dna_mutect2_vcf_extension="${source_dna_mutect2_vcf_extension:-${dna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
  source_rna_mutect2_vcf_extension="$(resolve_patient_placeholder "$source_rna_mutect2_vcf_extension" "$name")"
  source_dna_mutect2_vcf_extension="$(resolve_patient_placeholder "$source_dna_mutect2_vcf_extension" "$name")"
  source_rna_mutect2_vcf_extension="${source_rna_mutect2_vcf_extension%\}}"
  source_dna_mutect2_vcf_extension="${source_dna_mutect2_vcf_extension%\}}"

  rna_vcf_pref="${rna_vcf_dir}/${name}_${source_rna_mutect2_vcf_extension}"
  dna_vcf_pref="${dna_vcf_dir}/${name}_${source_dna_mutect2_vcf_extension}"
  rna_vcf="$(pick_first_existing \
    "$rna_vcf_pref" \
    "${rna_vcf_dir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${name}_${rna_tumor_label:-RNA_TUMOR}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${name}_RNA_TUMOR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${name}_RNA_TUMOUR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
  )" || rna_vcf=""
  dna_vcf="$(pick_first_existing \
    "$dna_vcf_pref" \
    "${dna_vcf_dir}/${name}_${dna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${name}_${dna_tumor_label:-DNA_TUMOR}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${name}_DNA_TUMOR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${name}_DNA_TUMOUR_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz" \
  )" || dna_vcf="$dna_vcf_pref"

  germ_merge_ext="${germline_for_merge_extension:-${gdna4_vcf_extension:-${output_extension_30:-3.0.Filtered.vcf}}}"
  germ_vcf="${germ_dir}/${name}_${germ_merge_ext}"

  dna_bam="${bamdir}/${name}_${dna_label}/${name}_${dna_bam_suffix:-md.bam}"

  outdir="${dna_vcf_dir}"
  dna_only_merged_vcf_extension="${dna_only_merged_vcf_extension:-dna_only6.DNAt_DNAn_merged.vcf.gz}"
  dna_only_genotyped_vcf_extension="${dna_only_genotyped_vcf_extension:-dna_only7.DNAt_DNAn_merged_genotyped.vcf.gz}"
  dna_only_phased_vcf_extension="${dna_only_phased_vcf_extension:-dna_only7.DNAt_DNAn_merged_phased.vcf.gz}"

  merged_vcf="${outdir}/${name}_${dna_only_merged_vcf_extension}"
  genotyped_vcf="${outdir}/${name}_${dna_only_genotyped_vcf_extension}"
  phased_vcf="${outdir}/${name}_${dna_only_phased_vcf_extension}"
  run_log="${phased_vcf%.vcf.gz}.dna_only.log.txt"

  if [ "$allow_with_rna" -eq 0 ] && [ -n "$rna_vcf" ]; then
    echo "[skip] ${prefix}.${name}: RNA source VCF exists, run RNA pipeline instead: $rna_vcf"
    continue
  fi

  if [ "$force" -eq 0 ] && [ -f "$phased_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $phased_vcf (use -f to overwrite)"
    continue
  fi

  if ! precheck_vcf_or_vcfgz "$germ_vcf" "${prefix}.${name}"; then
    continue
  fi
  if ! precheck_vcf_or_vcfgz "$dna_vcf" "${prefix}.${name}"; then
    continue
  fi
  if [ ! -f "$dna_bam" ]; then
    echo "ERROR: missing DNA BAM: $dna_bam" >&2
    echo "[skip] ${prefix}.${name}: missing DNA BAM" >&2
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
SCRIPT
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'germ_vcf=%q\n' "$germ_vcf"
    printf 'dna_vcf=%q\n' "$dna_vcf"
    printf 'dna_bam=%q\n' "$dna_bam"
    printf 'FASTA=%q\n' "$FASTA"
    printf 'DICT=%q\n' "$DICT"
    printf 'merged_vcf=%q\n' "$merged_vcf"
    printf 'genotyped_vcf=%q\n' "$genotyped_vcf"
    printf 'phased_vcf=%q\n' "$phased_vcf"
    printf 'out_normal=%q\n' "$out_normal_label"
    printf 'out_dna=%q\n' "$dna_label"
    printf 'run_log=%q\n' "$run_log"
    cat <<'SCRIPT'
mkdir -p "$(dirname "$merged_vcf")"

bash "$rnae_scripts/rnae6_merge_germ_som_dna_vcf.sh" \
  -g "$germ_vcf" -d "$dna_vcf" -o "$merged_vcf" --out-normal "$out_normal" --out-dna "$out_dna" 2>&1 | tee "$run_log"

bash "$rnae_scripts/rnae7_genotype_and_phase_dna_merged_vcf.sh" \
  --merged "$merged_vcf" --dna-bam "$dna_bam" --fasta "$FASTA" --dict "$DICT" \
  --out-genotyped "$genotyped_vcf" --out-phased "$phased_vcf" \
  --normal-name "$out_normal" --dna-name "$out_dna" 2>&1 | tee -a "$run_log"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=8,mem=48gb,walltime="01:00:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"
  echo ".. logs and reports saved in $scriptdir"
  sleep 0.3
done < "$samples"
