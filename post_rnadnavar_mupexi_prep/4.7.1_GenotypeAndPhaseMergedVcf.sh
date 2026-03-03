#!/usr/bin/bash
set -euo pipefail

# 4.7.1 Genotype union alleles back into BAMs and phase merged VCF.

usage() {
  cat <<'USAGE'
Usage: bash 4.7.1_GenotypeAndPhaseMergedVcf.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${bamdir:?CONFIG must define bamdir}"
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${FASTA:?CONFIG must define FASTA}"
: "${DICT:?CONFIG must define DICT}"
: "${merged_vcf_extension:?CONFIG must define merged_vcf_extension}"
: "${genotyped_vcf_extension:?CONFIG must define genotyped_vcf_extension}"
: "${phased_vcf_extension:?CONFIG must define phased_vcf_extension}"
: "${rna_bam_smfixed_suffix:?CONFIG must define rna_bam_smfixed_suffix}"

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
  echo "Running 4.7.1 for all samples in $samples"
else
  echo "Running 4.7.1 only for $sample"
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

  merged_vcf="${outdir_only}/${name}_${merged_vcf_extension}"
  genotyped_vcf="${outdir_only}/${name}_${genotyped_vcf_extension}"
  phased_vcf="${outdir_only}/${name}_${phased_vcf_extension}"

  dna_dir_label="${dna_tumor_label:-DNA_TUMOR}"
  rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"

  dna_bam="${bamdir}/${name}_${dna_dir_label}/${name}_${dna_dir_label}.md.bam"
  rna_bam_fixed="${bamdir}/${name}_${rna_dir_label}/${name}_${rna_bam_smfixed_suffix}"
  rna_bam_default="${bamdir}/${name}_${rna_dir_label}/${name}_${rna_dir_label}.md.bam"
  if [ -f "$rna_bam_fixed" ]; then
    rna_bam="$rna_bam_fixed"
  else
    rna_bam="$rna_bam_default"
  fi

  run_log="${phased_vcf%.vcf.gz}.genotype_and_phase.log.txt"

  if [ "$force" -eq 0 ] && [ -f "$phased_vcf" ]; then
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load ngs tools htslib/1.23 bcftools/1.23 java/17-openjdk gatk/4.4.0.0 anaconda3/2025.06-1
SCRIPT
    printf 'merged_vcf=%q\n' "$merged_vcf"
    printf 'dna_bam=%q\n' "$dna_bam"
    printf 'rna_bam=%q\n' "$rna_bam"
    printf 'fasta=%q\n' "$FASTA"
    printf 'dict=%q\n' "$DICT"
    printf 'genotyped_vcf=%q\n' "$genotyped_vcf"
    printf 'phased_vcf=%q\n' "$phased_vcf"
    printf 'run_log=%q\n' "$run_log"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'out_dna_normal_label=%q\n' "${out_dna_normal_label:-DNA_NORMAL}"
    printf 'out_dna_tumor_label=%q\n' "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
    printf 'out_rna_tumor_label=%q\n' "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    cat <<'SCRIPT'
if [ ! -f "$merged_vcf" ]; then
  echo "ERROR: missing merged VCF (run 4.6 first): $merged_vcf" >&2
  exit 1
fi
if [ ! -f "$dna_bam" ]; then
  echo "ERROR: missing DNA tumour BAM: $dna_bam" >&2
  exit 1
fi
if [ ! -f "$rna_bam" ]; then
  echo "ERROR: missing RNA tumour BAM: $rna_bam" >&2
  exit 1
fi

bash "$rnae_scripts/rnae7_genotype_and_phase_merged_vcf.sh" \
  --merged "$merged_vcf" \
  --dna-bam "$dna_bam" \
  --rna-bam "$rna_bam" \
  --fasta "$fasta" \
  --dict "$dict" \
  --out-genotyped "$genotyped_vcf" \
  --out-phased "$phased_vcf" \
  --threads 8 \
  --normal-name "$out_dna_normal_label" \
  --dna-name "$out_dna_tumor_label" \
  --rna-name "$out_rna_tumor_label" 2>&1 | tee "$run_log"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=8,mem=48gb,walltime="01:00:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
