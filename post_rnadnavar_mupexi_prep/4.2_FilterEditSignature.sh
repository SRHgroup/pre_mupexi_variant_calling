#!/usr/bin/bash
set -euo pipefail

# 4.2 Split multi-allelics and filter to RNA-editing signatures (ADAR/APOBEC3).

usage() {
  cat <<'USAGE'
Usage: bash 4.2_FilterEditSignature.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${rna_only_vcf_extension:?CONFIG must define rna_only_vcf_extension}"
: "${split_vcf_extension:?CONFIG must define split_vcf_extension}"
: "${filtered_edit_labeled_vcf_extension:?CONFIG must define filtered_edit_labeled_vcf_extension}"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}"
    "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do value="${value%_${label}}"; done
  printf '%s\n' "$value"
}

if [ -z "${sample:-}" ]; then
  echo "Running 4.2 for all samples in $samples"
else
  echo "Running 4.2 only for $sample"
fi

prefix=$(basename "${BASH_SOURCE[0]}" .sh)
scriptdir="${vcfdir}/${prefix}.logs_and_reports"
logdir="${scriptdir}/logs"
repdir="${scriptdir}/reports"
mkdir -p "$logdir" "$repdir"

while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,	 ]+' '{print $1}')
  [ -z "${sample:-}" ] || [ "$sample_name" = "$sample" ] || continue

  name=$(sample_base_name "$sample_name")
  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  rna_only_vcf="${outdir_only}/${name}_${rna_only_vcf_extension}"
  split_vcf="${outdir_only}/${name}_${split_vcf_extension}"
  filtered_vcf="${outdir_only}/${name}_${filtered_edit_labeled_vcf_extension}"
  stats_file="${filtered_vcf%.vcf.gz}.stats.tsv"
  counts_file="${filtered_vcf%.vcf.gz}.counts.txt"

  if [ "$force" -eq 0 ] && [ -f "$filtered_vcf" ]; then
    continue
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
  echo "ERROR: missing input (run 4.1 first): $rna_only_vcf" >&2
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

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=4,mem=12gb,walltime="00:10:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
