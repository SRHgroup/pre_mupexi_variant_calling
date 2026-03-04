#!/usr/bin/bash
set -euo pipefail

# rna2: split multi-allelics and filter to RNA-editing signatures (ADAR/APOBEC3).

usage() {
  cat <<'USAGE'
Usage: bash rna2_FilterEditSignature.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${vcfdir:?CONFIG must define vcfdir}"
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${rna1_vcf_extension:?CONFIG must define rna1_vcf_extension}"
: "${rna2_split_vcf_extension:?CONFIG must define rna2_split_vcf_extension}"
: "${rna2_labeled_vcf_extension:?CONFIG must define rna2_labeled_vcf_extension}"

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

precheck_vcfgz() {
  local path="$1"
  local tag="$2"
  local hint="$3"
  if [ ! -f "$path" ]; then
    echo "[precheck] ${tag}: missing input VCF: $path (${hint})" >&2
    return 1
  fi
  if ! gzip -t "$path" 2>/dev/null; then
    echo "[precheck] ${tag}: invalid gzip VCF: $path" >&2
    return 1
  fi
  # Primary fast check
  if zgrep -m1 '^#CHROM[[:space:]]' "$path" >/dev/null 2>&1; then
    return 0
  fi
  # Fallback: inspect first non-meta line (handles odd grep behavior/BOM/CRLF)
  local first_non_meta
  first_non_meta="$(gzip -dc "$path" 2>/dev/null | awk '
    BEGIN{found=0}
    /^##/ {next}
    {
      sub(/\r$/,"",$0)
      sub(/^\xef\xbb\xbf/,"",$0)
      print
      found=1
      exit
    }
    END{ if (!found) exit 2 }
  ' || true)"
  if [[ "$first_non_meta" =~ ^#CHROM[[:space:]] ]]; then
    echo "[precheck] ${tag}: header check recovered via fallback parser: $path" >&2
    return 0
  fi
  if ! zgrep -m1 '^#CHROM' "$path" >/dev/null 2>&1; then
    echo "[precheck] ${tag}: malformed VCF header (missing #CHROM): $path" >&2
    return 1
  fi
}

if [ -z "${sample:-}" ]; then
  echo "Running rna2 for all samples in $samples"
else
  echo "Running rna2 only for $sample"
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
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,	 ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1
  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  rna_only_vcf="${outdir_only}/${name}_${rna1_vcf_extension}"
  split_vcf="${outdir_only}/${name}_${rna2_split_vcf_extension}"
  filtered_vcf="${outdir_only}/${name}_${rna2_labeled_vcf_extension}"
  stats_file="${filtered_vcf%.vcf.gz}.stats.tsv"
  counts_file="${filtered_vcf%.vcf.gz}.counts.txt"

  if [ "$force" -eq 0 ] && [ -f "$filtered_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $filtered_vcf (use -f to overwrite)"
    continue
  fi
  if ! precheck_vcfgz "$rna_only_vcf" "${prefix}.${name}" "run rna1 first"; then
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
module load ngs tools htslib/1.23 bcftools/1.23 gatk/4.4.0.0 anaconda3/2025.06-1
SCRIPT
    printf 'outdir_only=%q\n' "$outdir_only"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'rna_only_vcf=%q\n' "$rna_only_vcf"
    printf 'split_vcf=%q\n' "$split_vcf"
    printf 'filtered_vcf=%q\n' "$filtered_vcf"
    printf 'stats_file=%q\n' "$stats_file"
    printf 'counts_file=%q\n' "$counts_file"
    cat <<'SCRIPT'
mkdir -p "$outdir_only"

if [ ! -f "$rna_only_vcf" ]; then
  echo "ERROR: missing input (run rna1 first): $rna_only_vcf" >&2
  exit 1
fi

bcftools norm -m -both -Oz -o "$split_vcf" "$rna_only_vcf"
bcftools index -t "$split_vcf"

python3 "$rnae_scripts/rnae2_annotate_and_filter_edit_sig.py" -i "$split_vcf" -o "$filtered_vcf" 2> "$stats_file"
bcftools index -t "$filtered_vcf"

{
  echo "Counts in filtered VCF"
  echo -n "total kept:     "; bcftools view -H "$filtered_vcf" | wc -l
  echo -n "ADAR_SIG:       "; bcftools view -H -i 'INFO/ADAR_SIG' "$filtered_vcf" | wc -l
  echo -n "APOBEC3_SIG:    "; bcftools view -H -i 'INFO/APOBEC3_SIG' "$filtered_vcf" | wc -l
} > "$counts_file"

SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=4,mem=12gb,walltime="00:10:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
