#!/usr/bin/bash
set -euo pipefail

# 4.6 Merge three VCFs: DNA normal germline + DNA tumour somatic + RNA tumour editing.

usage() {
  cat <<'USAGE'
Usage: bash 4.6_MergeDnaRnaVcfs.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${rna_vcf_knownsites_extension:?CONFIG must define rna_vcf_knownsites_extension}"
: "${merged_vcf_extension:?CONFIG must define merged_vcf_extension}"

# Germline VCF suffix to merge. Defaults to 3.0 output if not set.
ndna_vcf="${ndna_vcf:-${output_extension_30:-3.0.Filtered.vcf}}"

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

if [ -z "${sample:-}" ]; then
  echo "Running 4.6 for all samples in $samples"
else
  echo "Running 4.6 only for $sample"
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
  dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
  germline_dir="${vcfdir}/${name}_${out_normal_label}"

  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"
  dna_dir="${vcfdir}/${name}_${dna_label}_vs_${name}_${out_normal_label}"
  germ_in="${germline_dir}/${name}_${ndna_vcf}"
  dna_tumour_vcf="${dna_dir}/${name}_${dna_label}_vs_${name}_${out_normal_label}.mutect2.filtered.vcf.gz"
  rna_vcf="${outdir_only}/${name}_${rna_vcf_knownsites_extension}"

  merged_vcf="${outdir_only}/${name}_${merged_vcf_extension}"
  merge_log="${merged_vcf%.vcf.gz}.merge.log.txt"

  if [ "$force" -eq 0 ] && [ -f "$merged_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $merged_vcf (use -f to overwrite)"
    continue
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
module load ngs tools htslib/1.23 bcftools/1.23 anaconda3/2025.06-1
SCRIPT
    printf 'outdir_only=%q\n' "$outdir_only"
    printf 'dna_dir=%q\n' "$dna_dir"
    printf 'germline_dir=%q\n' "$germline_dir"
    printf 'name=%q\n' "$name"
    printf 'ndna_vcf=%q\n' "$ndna_vcf"
    printf 'germ_in=%q\n' "$germ_in"
    printf 'dna_tumour_vcf=%q\n' "$dna_tumour_vcf"
    printf 'rna_vcf=%q\n' "$rna_vcf"
    printf 'merged_vcf=%q\n' "$merged_vcf"
    printf 'merge_log=%q\n' "$merge_log"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    printf 'out_dna_normal_label=%q\n' "${out_dna_normal_label:-DNA_NORMAL}"
    printf 'out_dna_tumor_label=%q\n' "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
    printf 'out_rna_tumor_label=%q\n' "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    cat <<'SCRIPT'
if [ ! -f "$dna_tumour_vcf" ]; then
  echo "ERROR: missing DNA tumour VCF: $dna_tumour_vcf" >&2
  exit 1
fi
if [ ! -f "$rna_vcf" ]; then
  echo "ERROR: missing RNA VCF (run 4.5 first): $rna_vcf" >&2
  exit 1
fi

if [ -f "$germ_in.gz" ]; then
  germ_vcf="$germ_in.gz"
elif echo "$germ_in" | grep -q '\.vcf\.gz$'; then
  germ_vcf="$germ_in"
elif [ -f "$germ_in" ]; then
  germ_vcf="$germline_dir/${name}_${ndna_vcf}.gz"
  bgzip -c "$germ_in" > "$germ_vcf"
else
  echo "ERROR: missing germline VCF: $germ_in (or $germ_in.gz)" >&2
  exit 1
fi

bcftools index -f -t "$dna_tumour_vcf" || true
bcftools index -f -t "$rna_vcf" || true
bcftools index -f -t "$germ_vcf" || true

bash "$rnae_scripts/rnae6_merge_germ_som_rna_vcf.sh" \
  -g "$germ_vcf" -d "$dna_tumour_vcf" -r "$rna_vcf" -o "$merged_vcf" --rna-keep-all \
  --out-normal "$out_dna_normal_label" --out-dna "$out_dna_tumor_label" --out-rna "$out_rna_tumor_label" 2>&1 | tee "$merge_log"

bcftools index -t "$merged_vcf"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=4,mem=24gb,walltime="00:16:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
