#!/usr/bin/bash
set -euo pipefail

# 4.5 Filter summarised RNA VCF by AF/DP/ALT-reads thresholds.

usage() {
  cat <<'USAGE'
Usage: bash 4.5_FilterByAfDpAr.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${rna_summarised_vcf_extension:?CONFIG must define rna_summarised_vcf_extension}"
: "${rna_vcf_knownsites_extension:?CONFIG must define rna_vcf_knownsites_extension}"

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

precheck_input_vcf() {
  local path="$1"
  local tag="$2"
  if [ ! -f "$path" ]; then
    echo "[precheck] ${tag}: missing input VCF (run 4.4 first): $path" >&2
    return 1
  fi
  if ! gzip -t "$path" 2>/dev/null; then
    echo "[precheck] ${tag}: input VCF is not valid gzip: $path" >&2
    return 1
  fi
  if ! zgrep -m1 '^#CHROM' "$path" >/dev/null 2>&1; then
    echo "[precheck] ${tag}: input VCF header missing #CHROM (corrupt/malformed): $path" >&2
    return 1
  fi
  return 0
}

if [ -z "${sample:-}" ]; then
  echo "Running 4.5 for all samples in $samples"
else
  echo "Running 4.5 only for $sample"
fi

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
  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  in_vcf="${outdir_only}/${name}_${rna_summarised_vcf_extension}"
  out_vcf="${outdir_only}/${name}_${rna_vcf_knownsites_extension}"
  stats_file="${out_vcf%.vcf.gz}.filter_stats.tsv"
  counts_file="${out_vcf%.vcf.gz}.counts.txt"

  if [ "$force" -eq 0 ] && [ -f "$out_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $out_vcf (use -f to overwrite)"
    continue
  fi
  if ! precheck_input_vcf "$in_vcf" "${prefix}.${name}"; then
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load ngs tools htslib/1.23 bcftools/1.23 gatk/4.4.0.0 anaconda3/2025.06-1
SCRIPT
    printf 'in_vcf=%q\n' "$in_vcf"
    printf 'out_vcf=%q\n' "$out_vcf"
    printf 'stats_file=%q\n' "$stats_file"
    printf 'counts_file=%q\n' "$counts_file"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'dna_normal_label=%q\n' "${dna_normal_label:-DNA_NORMAL}"
    printf 'rna_tumor_label=%q\n' "${rna_tumor_label:-RNA_TUMOR}"
    printf 'out_dna_normal_label=%q\n' "${out_dna_normal_label:-DNA_NORMAL}"
    printf 'out_rna_tumor_label=%q\n' "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    cat <<'SCRIPT'
if [ ! -f "$in_vcf" ]; then
  echo "ERROR: missing input (run 4.4 first): $in_vcf" >&2
  exit 1
fi

if ! gzip -t "$in_vcf" 2>/dev/null; then
  echo "ERROR: 4.5 input is not a valid gzip VCF: $in_vcf" >&2
  echo "Re-run step 4.4 with -f for this sample to regenerate the file." >&2
  exit 1
fi
if ! bcftools view -h "$in_vcf" | grep -q '^#CHROM'; then
  echo "ERROR: 4.5 input header is malformed (missing #CHROM): $in_vcf" >&2
  echo "Re-run step 4.4 with -f for this sample to regenerate the file." >&2
  exit 1
fi

normal_sample=$(bcftools query -l "$in_vcf" | grep -x "$out_dna_normal_label" | head -n 1 || true)
rna_sample=$(bcftools query -l "$in_vcf" | grep -x "$out_rna_tumor_label" | head -n 1 || true)

if [ -z "$normal_sample" ]; then
  normal_sample=$(bcftools query -l "$in_vcf" | grep -E "(${dna_normal_label})$" | head -n 1 || true)
fi
if [ -z "$rna_sample" ]; then
  a="$rna_tumor_label"
  alt="$a"
  if echo "$a" | grep -q TUMOR; then alt="${a/TUMOR/TUMOUR}"; fi
  if echo "$a" | grep -q TUMOUR; then alt="${a/TUMOUR/TUMOR}"; fi
  rna_sample=$(bcftools query -l "$in_vcf" | grep -E "(_)?(${a}|${alt})$" | head -n 1 || true)
  if [ -z "$rna_sample" ]; then
    rna_sample=$(bcftools query -l "$in_vcf" | grep -E "(_)?(${a}|${alt})(\\.|$)" | head -n 1 || true)
  fi
fi

if [ -z "$normal_sample" ] || [ -z "$rna_sample" ]; then
  echo "ERROR: could not detect required samples in $in_vcf" >&2
  echo "Detected normal='$normal_sample' rna='$rna_sample'" >&2
  echo "Header samples:" >&2
  bcftools query -l "$in_vcf" >&2 || true
  exit 1
fi

python3 "$rnae_scripts/rnae5_filter_by_AF_DP_AR.py" \
  -i "$in_vcf" -o "$out_vcf" \
  --normal-sample "$normal_sample" --rna-sample "$rna_sample" \
  --stats 2> "$stats_file"

if [ ! -s "$out_vcf" ]; then
  echo "ERROR: 4.5 produced empty output file: $out_vcf" >&2
  exit 1
fi
if ! bcftools view -h "$out_vcf" | grep -q '^#CHROM'; then
  echo "ERROR: 4.5 produced malformed VCF (missing #CHROM): $out_vcf" >&2
  rm -f "$out_vcf" "$out_vcf.tbi"
  exit 1
fi

bcftools index -t "$out_vcf"

{
  echo "Counts in filtered VCF"
  echo -n "input_total:  "; bcftools view -H "$in_vcf" | wc -l
  echo -n "output_total: "; bcftools view -H "$out_vcf" | wc -l
} > "$counts_file"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=4,mem=16gb,walltime="00:08:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
