#!/usr/bin/bash
set -euo pipefail

# gdna4: filter selected germline calls by adjacency to somatic variants (optional RNA variants).

usage() {
  cat <<'USAGE'
Usage: bash gdna4_FilterGermlineByAdjacency.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
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
: "${gdna3_vcf_extension:?CONFIG must define gdna3_vcf_extension}"
: "${gdna4_vcf_extension:?CONFIG must define gdna4_vcf_extension}"

distance="${germline_proximity_distance:-33}"
# RNA proximity behavior:
#   auto (default): include RNA proximity if an RNA VCF is found
#   1/true/yes/on: require attempt, include when found, warn when missing
#   0/false/no/off: disable RNA proximity
include_rna_mode="${germline_include_rna_proximity:-auto}"

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

helpers_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"

detect_rna_vcf_for_patient() {
  local patient="$1"
  local rna_dir="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}"
  local path=""

  # Optional explicit absolute/relative path template:
  #   germline_rna_proximity_vcf="/path/to/.../{patient}_4.5....vcf.gz"
  if [ -n "${germline_rna_proximity_vcf:-}" ]; then
    path="${germline_rna_proximity_vcf//\{patient\}/$patient}"
    if [ -f "$path" ]; then
      printf '%s\n' "$path"
      return 0
    fi
  fi

  # Optional extension override for files inside the patient RNA_vs_normal folder.
  if [ -n "${germline_rna_proximity_vcf_extension:-}" ]; then
    path="${rna_dir}/${patient}_${germline_rna_proximity_vcf_extension}"
    if [ -f "$path" ]; then
      printf '%s\n' "$path"
      return 0
    fi
  fi

  # Default candidate order: prefer filtered RNA-editing sets, then source mutect2 RNA callset.
  for ext in \
    "${rna5_qced_vcf_extension:-}" \
    "${rna4_summarised_vcf_extension:-}" \
    "${rna3_knownsites_vcf_extension:-}" \
    "${rna2_labeled_vcf_extension:-}" \
    "${rna1_vcf_extension:-}" \
    "${out_rna_label}_vs_${out_normal_label}.mutect2.filtered.vcf.gz"; do
    [ -n "$ext" ] || continue
    if [[ "$ext" = *.vcf.gz || "$ext" = *.vcf ]]; then
      path="${rna_dir}/${patient}_${ext}"
    else
      path="${rna_dir}/${patient}_${ext}.vcf.gz"
    fi
    if [ -f "$path" ]; then
      printf '%s\n' "$path"
      return 0
    fi
  done

  printf '\n'
  return 0
}

rna_mode_enabled() {
  case "${include_rna_mode,,}" in
    1|true|yes|on|auto) return 0 ;;
    0|false|no|off) return 1 ;;
    *)
      echo "ERROR: invalid germline_include_rna_proximity='${include_rna_mode}' (use auto|1|0)" >&2
      exit 1
      ;;
  esac
}

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
  germ_in="${germline_dir}/${name}_${gdna3_vcf_extension}"
  germ_out="${germline_dir}/${name}_${gdna4_vcf_extension}"

  dna_dir="${vcfdir}/${name}_${out_dna_label}_vs_${name}_${out_normal_label}"
  somatic_vcf="${dna_dir}/${name}_${out_dna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz"

  rna_vcf=""
  if rna_mode_enabled; then
    rna_vcf="$(detect_rna_vcf_for_patient "$name")"
    if [ -n "$rna_vcf" ]; then
      echo "[info] ${prefix}.${name}: RNA proximity enabled with $rna_vcf"
    elif [[ "${include_rna_mode,,}" = "auto" ]]; then
      echo "[info] ${prefix}.${name}: RNA proximity auto-mode found no RNA VCF; using DNA-only proximity"
    else
      echo "[warn] ${prefix}.${name}: RNA proximity requested but no RNA VCF found; using DNA-only proximity"
    fi
  else
    echo "[info] ${prefix}.${name}: RNA proximity disabled by config"
  fi

  if [ "$force" -eq 0 ] && [ -f "$germ_out" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $germ_out (use -f to overwrite)"
    continue
  fi
  if [ ! -s "$germ_in" ]; then
    echo "[precheck] ${prefix}.${name}: missing/empty germline input from gdna3: $germ_in" >&2
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi
  if [ ! -s "$somatic_vcf" ]; then
    echo "[precheck] ${prefix}.${name}: missing/empty DNA somatic VCF: $somatic_vcf" >&2
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
  echo "ERROR: missing germline VCF from gdna3: $germ_in" >&2
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

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=4,mem=8gb,walltime="00:04:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
