#!/usr/bin/bash
set -euo pipefail

if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  source "$PIPELINE_DEFAULTS"
fi
module load ${modules_dna_only_phase:-ngs tools htslib/1.23 bcftools/1.23 gatk/4.5.0.0 anaconda3/2025.06-1}

threads=8
merged_vcf=""
dna_bam=""
fasta=""
dict=""
out_genotyped=""
out_phased=""

normal_name="DNA_NORMAL"
dna_name="DNA_TUMOR"

usage() {
  cat <<USAGE
Usage:
  $0 --merged merged.vcf.gz --dna-bam DNA_TUMOR.bam \
     --fasta ref.fa.gz --dict ref.dict \
     --out-genotyped out.genotyped.vcf.gz --out-phased out.phased.vcf.gz [--threads 8]

Optional sample names:
  --normal-name  (default: DNA_NORMAL)
  --dna-name     (default: DNA_TUMOR)
USAGE
}

die(){ echo "[error] $*" >&2; exit 1; }

resolve_sample_name () {
  local requested="$1"
  shift
  local -a list=("$@")
  local s alt
  for s in "${list[@]}"; do
    [[ "$s" == "$requested" ]] && { echo "$s"; return 0; }
  done
  for s in "${list[@]}"; do
    [[ "${s^^}" == "${requested^^}" || "${s^^}" =~ _${requested^^}$ ]] && { echo "$s"; return 0; }
  done
  if [[ "$requested" == *TUMOR* ]]; then alt="${requested/TUMOR/TUMOUR}"; else alt="${requested/TUMOUR/TUMOR}"; fi
  for s in "${list[@]}"; do
    [[ "${s^^}" == "${alt^^}" || "${s^^}" =~ _${alt^^}$ || "${s^^}" =~ _${alt^^}(\.[0-9]+)?$ ]] && { echo "$s"; return 0; }
  done
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged) merged_vcf="$2"; shift 2;;
    --dna-bam) dna_bam="$2"; shift 2;;
    --fasta) fasta="$2"; shift 2;;
    --dict) dict="$2"; shift 2;;
    --out-genotyped) out_genotyped="$2"; shift 2;;
    --out-phased) out_phased="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
    --normal-name) normal_name="$2"; shift 2;;
    --dna-name) dna_name="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -n "$merged_vcf" ]] || die "Missing --merged"
[[ -n "$dna_bam" ]] || die "Missing --dna-bam"
[[ -n "$fasta" ]] || die "Missing --fasta"
[[ -n "$dict" ]] || die "Missing --dict"
[[ -n "$out_genotyped" ]] || die "Missing --out-genotyped"
[[ -n "$out_phased" ]] || die "Missing --out-phased"

[[ -f "$merged_vcf" ]] || die "Missing merged VCF: $merged_vcf"
[[ -f "$dna_bam" ]] || die "Missing DNA BAM: $dna_bam"
[[ -f "$fasta" ]] || die "Missing FASTA: $fasta"
[[ -f "${fasta}.fai" ]] || die "Missing FASTA index: ${fasta}.fai"
[[ -f "$dict" ]] || die "Missing dict: $dict"

command -v whatshap >/dev/null 2>&1 || die "whatshap not found"
command -v gatk >/dev/null 2>&1 || die "gatk not found"
command -v bcftools >/dev/null 2>&1 || die "bcftools not found"
command -v bgzip >/dev/null 2>&1 || die "bgzip not found"

bcftools index -f -t "$merged_vcf" >/dev/null 2>&1 || true
mapfile -t samples_arr < <(bcftools query -l "$merged_vcf")
[[ ${#samples_arr[@]} -ge 2 ]] || die "Merged VCF must contain at least normal + tumor samples"

normal_in=$(resolve_sample_name "$normal_name" "${samples_arr[@]}") || die "$normal_name not found in merged VCF"
dna_in=$(resolve_sample_name "$dna_name" "${samples_arr[@]}") || die "$dna_name not found in merged VCF"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

union_sites_vcf="${tmpdir}/union.sites.vcf.gz"
bcftools view -G -Oz -o "$union_sites_vcf" "$merged_vcf"
bcftools index -f -t "$union_sites_vcf"

intervals_bed="${tmpdir}/union.sites.bed"
bcftools query -f '%CHROM\t%POS\n' "$union_sites_vcf" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' > "$intervals_bed"

dna_raw="${tmpdir}/dna.raw.vcf.gz"
gatk --java-options "-Xmx16g" HaplotypeCaller \
  -R "$fasta" -I "$dna_bam" -L "$intervals_bed" \
  --alleles "$union_sites_vcf" \
  --output-mode EMIT_ALL_ACTIVE_SITES \
  --disable-tool-default-annotations \
  --active-probability-threshold 0.0 \
  --annotations-to-exclude TandemRepeat \
  --native-pair-hmm-threads "$threads" \
  -O "$dna_raw"
bcftools index -f -t "$dna_raw" >/dev/null 2>&1 || true

dna_fix="${tmpdir}/${dna_name}.vcf.gz"
echo -e "$(bcftools query -l "$dna_raw" | head -n1)\t${dna_name}" > "${tmpdir}/dna.map"
bcftools reheader -s "${tmpdir}/dna.map" -o "$dna_fix" "$dna_raw"
bcftools index -f -t "$dna_fix"

normal_vcf="${tmpdir}/${normal_name}.vcf.gz"
bcftools view -s "$normal_in" -Oz -o "$normal_vcf" "$merged_vcf"
bcftools index -f -t "$normal_vcf"

info_rules="SOURCE_SET:join"
bcftools merge --threads "$threads" -m all --info-rules "$info_rules" -Oz -o "$out_genotyped" "$normal_vcf" "$dna_fix"
bcftools index -f -t "$out_genotyped"

echo "[info] Wrote genotyped DNA-only VCF: $out_genotyped"

whatshap phase --reference "$fasta" --sample "$dna_name" --ignore-read-groups \
  --output /dev/stdout "$out_genotyped" "$dna_bam" | bgzip -c > "$out_phased"
bcftools index -f -t "$out_phased"

echo "[info] Wrote phased DNA-only VCF: $out_phased"
