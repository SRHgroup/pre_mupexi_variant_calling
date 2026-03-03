#!/usr/bin/bash
set -euo pipefail

# 4.3 Annotate filtered RNA VCF with known RNA-editing sites database (INFO/KNOWN_RNAEDIT_DB).

usage() {
  cat <<'USAGE'
Usage: bash 4.3_AnnotateKnownSites.sh -c CONFIG [-s SAMPLE] [-f]
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
: "${knownsites:?CONFIG must define knownsites}"
: "${filtered_edit_labeled_vcf_extension:?CONFIG must define filtered_edit_labeled_vcf_extension}"
: "${annot_vcf_extension:?CONFIG must define annot_vcf_extension}"

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
  echo "Running 4.3 for all samples in $samples"
else
  echo "Running 4.3 only for $sample"
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

  in_vcf="${outdir_only}/${name}_${filtered_edit_labeled_vcf_extension}"
  annot_vcf="${outdir_only}/${name}_${annot_vcf_extension}"
  annot_log="${annot_vcf%.vcf.gz}.counts.tsv"

  if [ "$force" -eq 0 ] && [ -f "$annot_vcf" ]; then
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
    printf 'in_vcf=%q\n' "$in_vcf"
    printf 'annot_vcf=%q\n' "$annot_vcf"
    printf 'annot_log=%q\n' "$annot_log"
    printf 'knownsites=%q\n' "$knownsites"
    cat <<'SCRIPT'
mkdir -p "$outdir_only"

if [ ! -f "$in_vcf" ]; then
  echo "ERROR: missing input (run 4.2 first): $in_vcf" >&2
  exit 1
fi
if [ ! -f "$knownsites" ]; then
  echo "ERROR: missing knownsites DB: $knownsites" >&2
  exit 1
fi

hdr_tmp=$(mktemp)
trap 'rm -f "$hdr_tmp"' EXIT
cat > "$hdr_tmp" <<'HDR'
##INFO=<ID=KNOWN_RNAEDIT_DB,Number=.,Type=String,Description="Known RNA-editing DB hit(s); merged sources keyed by CHROM,POS,REF,ALT">
HDR

bcftools annotate -a "$knownsites" -c CHROM,POS,REF,ALT,INFO/KNOWN_RNAEDIT_DB -h "$hdr_tmp" -Oz -o "$annot_vcf" "$in_vcf"
bcftools index -t "$annot_vcf"

in_total=$(bcftools view -H "$in_vcf" | wc -l)
out_total=$(bcftools view -H "$annot_vcf" | wc -l)
out_hits=$(bcftools view -H -i 'INFO/KNOWN_RNAEDIT_DB!=""' "$annot_vcf" | wc -l)

redi_hits=$(bcftools view -H -i 'INFO/KNOWN_RNAEDIT_DB~"REDIportal"' "$annot_vcf" | wc -l || true)
radar_hits=$(bcftools view -H -i 'INFO/KNOWN_RNAEDIT_DB~"RADAR"' "$annot_vcf" | wc -l || true)
asaoka_hits=$(bcftools view -H -i 'INFO/KNOWN_RNAEDIT_DB~"Asaoka"' "$annot_vcf" | wc -l || true)

{
  printf "metric\tvalue\n"
  printf "input_total_variants\t%s\n" "$in_total"
  printf "output_total_variants\t%s\n" "$out_total"
  printf "output_knownsite_hits\t%s\n" "$out_hits"
  printf "hits_REDIportal\t%s\n" "$redi_hits"
  printf "hits_RADAR\t%s\n" "$radar_hits"
  printf "hits_Asaoka_APOBEC\t%s\n" "$asaoka_hits"
  printf "input_vcf\t%s\n" "$in_vcf"
  printf "knownsites_db\t%s\n" "$knownsites"
  printf "output_vcf\t%s\n" "$annot_vcf"
} > "$annot_log"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    -l nodes=1:ppn=4,mem=16gb,walltime="00:08:00:00" -r y -N "${prefix}.${name}" -o "$repdir" -e "$repdir" "$runscript"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
