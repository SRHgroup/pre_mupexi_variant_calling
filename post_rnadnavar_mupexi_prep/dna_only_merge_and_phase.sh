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

resolve_bam_candidates() {
  local patient="$1"
  local tumor_label="$2"
  local suffix="${dna_bam_suffix:-md.bam}"
  local alt_label="$tumor_label"
  if [[ "$tumor_label" == *TUMOR* ]]; then
    alt_label="${tumor_label/TUMOR/TUMOUR}"
  elif [[ "$tumor_label" == *TUMOUR* ]]; then
    alt_label="${tumor_label/TUMOUR/TUMOR}"
  fi
  printf '%s\n' \
    "${bamdir}/${patient}_${tumor_label}/${patient}_${tumor_label}.${suffix}" \
    "${bamdir}/${patient}_${tumor_label}/${patient}_${suffix}" \
    "${bamdir}/${patient}_${tumor_label}/${patient}_${tumor_label}.md.bam" \
    "${bamdir}/${patient}_${tumor_label}/${patient}_md.bam" \
    "${bamdir}/${patient}_${alt_label}/${patient}_${alt_label}.${suffix}" \
    "${bamdir}/${patient}_${alt_label}/${patient}_${suffix}" \
    "${bamdir}/${patient}_${alt_label}/${patient}_${alt_label}.md.bam" \
    "${bamdir}/${patient}_${alt_label}/${patient}_md.bam"
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

  dna_bam="$(pick_first_existing $(resolve_bam_candidates "$name" "$dna_label"))" || dna_bam=""
  if [ -z "${dna_bam:-}" ]; then
    dna_bam="${bamdir}/${name}_${dna_label}/${name}_${dna_label}.${dna_bam_suffix:-md.bam}"
  fi

  outdir="${dna_vcf_dir}"
  dna_only_merged_vcf_extension="${dna_only_merged_vcf_extension:-dna_only6.DNAt_DNAn_merged.vcf.gz}"
  dna_only_genotyped_vcf_extension="${dna_only_genotyped_vcf_extension:-dna_only7.DNAt_DNAn_merged_genotyped.vcf.gz}"
  dna_only_phased_vcf_extension="${dna_only_phased_vcf_extension:-dna_only7.DNAt_DNAn_merged_phased.vcf.gz}"

  merged_vcf="${outdir}/${name}_${dna_only_merged_vcf_extension}"
  genotyped_vcf="${outdir}/${name}_${dna_only_genotyped_vcf_extension}"
  phased_vcf="${outdir}/${name}_${dna_only_phased_vcf_extension}"
  merge_log="${merged_vcf%.vcf.gz}.merge.log.txt"
  phase_log="${phased_vcf%.vcf.gz}.genotype_and_phase.log.txt"

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

  merge_job_name="${prefix}.merge.${name}"
  phase_job_name="${prefix}.phase.${name}"
  if [ "${SKIP_RUNNING:-0}" = "1" ]; then
    active_merge=""
    active_phase=""
    if command -v qselect >/dev/null 2>&1; then
      active_merge="$(qselect -u "${USER:-$(whoami)}" -N "$merge_job_name" 2>/dev/null | head -n1 || true)"
      active_phase="$(qselect -u "${USER:-$(whoami)}" -N "$phase_job_name" 2>/dev/null | head -n1 || true)"
    fi
    if [ -z "$active_merge" ] && command -v qstat >/dev/null 2>&1; then
      active_merge="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="$merge_job_name" '$4==n {print $1; exit}')"
    fi
    if [ -z "$active_phase" ] && command -v qstat >/dev/null 2>&1; then
      active_phase="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="$phase_job_name" '$4==n {print $1; exit}')"
    fi
    if [ -n "$active_merge" ] || [ -n "$active_phase" ]; then
      echo "[skip] ${phase_job_name}: matching merge/phase job already active: ${active_merge:-none} ${active_phase:-none}"
      continue
    fi
  fi
  submit_marker="${logdir}/submitted.${phase_job_name}.jobid"
  if [ -f "$submit_marker" ]; then
    prev_jobid="$(head -n1 "$submit_marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ] && command -v qstat >/dev/null 2>&1 && qstat "$prev_jobid" >/dev/null 2>&1; then
      echo "[skip] ${phase_job_name}: job already queued/running: ${prev_jobid}"
      continue
    fi
  fi

  merge_runscript="${logdir}/run.${name}.${prefix}.merge.sh"
  phase_runscript="${logdir}/run.${name}.${prefix}.phase.sh"
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
    printf 'out_normal=%q\n' "$out_normal_label"
    printf 'out_dna=%q\n' "$dna_label"
    printf 'merge_log=%q\n' "$merge_log"
    cat <<'SCRIPT'
mkdir -p "$(dirname "$merged_vcf")"

bash "$rnae_scripts/rnae6_merge_germ_som_dna_vcf.sh" \
  -g "$germ_vcf" -d "$dna_vcf" -o "$merged_vcf" --out-normal "$out_normal" --out-dna "$out_dna" 2>&1 | tee "$merge_log"
SCRIPT
  } > "$merge_runscript"
  chmod +x "$merge_runscript"

  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
SCRIPT
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'dna_bam=%q\n' "$dna_bam"
    printf 'FASTA=%q\n' "$FASTA"
    printf 'DICT=%q\n' "$DICT"
    printf 'merged_vcf=%q\n' "$merged_vcf"
    printf 'genotyped_vcf=%q\n' "$genotyped_vcf"
    printf 'phased_vcf=%q\n' "$phased_vcf"
    printf 'out_normal=%q\n' "$out_normal_label"
    printf 'out_dna=%q\n' "$dna_label"
    printf 'phase_log=%q\n' "$phase_log"
    cat <<'SCRIPT'
if [ ! -f "$merged_vcf" ]; then
  echo "ERROR: missing DNA-only merged VCF (merge step failed?): $merged_vcf" >&2
  exit 1
fi

bash "$rnae_scripts/rnae7_genotype_and_phase_dna_merged_vcf.sh" \
  --merged "$merged_vcf" --dna-bam "$dna_bam" --fasta "$FASTA" --dict "$DICT" \
  --out-genotyped "$genotyped_vcf" --out-phased "$phased_vcf" \
  --normal-name "$out_normal" --dna-name "$out_dna" 2>&1 | tee "$phase_log"
SCRIPT
  } > "$phase_runscript"
  chmod +x "$phase_runscript"

  merge_jobid=""
  if [ "$force" -eq 1 ] || [ ! -f "$merged_vcf" ]; then
    merge_jobid="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
      "${qsub_depend_arg[@]}" \
      -l nodes=1:ppn=4,mem=16gb,walltime="00:06:00:00" -r y -N "$merge_job_name" -o "$repdir" -e "$repdir" "$merge_runscript")"
    echo "$merge_jobid"
    echo "[submit] ${merge_job_name}: jobid=${merge_jobid}"
  else
    echo "[skip] ${merge_job_name}: merged output already exists: $merged_vcf"
  fi

  phase_depend_arg=()
  if [ -n "$merge_jobid" ]; then
    phase_depend_arg=(-W "depend=afterok:${merge_jobid}")
  fi
  phase_jobid="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${phase_depend_arg[@]}" \
    -l nodes=1:ppn=8,mem=48gb,walltime="01:00:00:00" -r y -N "$phase_job_name" -o "$repdir" -e "$repdir" "$phase_runscript")"
  echo "$phase_jobid"
  printf '%s\n' "$phase_jobid" > "$submit_marker"
  echo "[submit] ${phase_job_name}: jobid=${phase_jobid}"
  echo ".. logs and reports saved in $scriptdir"
  sleep 0.3
done < "$samples"
