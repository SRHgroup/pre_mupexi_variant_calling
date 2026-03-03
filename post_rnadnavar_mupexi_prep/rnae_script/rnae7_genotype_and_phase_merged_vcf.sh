#!/usr/bin/bash
set -euo pipefail

# rnae7: genotype union alleles (from merged VCF) back into BAMs and phase with whatshap.
# Sample-name aware (supports TUMOR/TUMOUR and configurable output names).

threads=8
merged_vcf=""
dna_bam=""
rna_bam=""
fasta=""
dict=""
out_genotyped=""
out_phased=""

normal_name="DNA_NORMAL"
dna_name="DNA_TUMOR"
rna_name="RNA_TUMOR"

usage() {
  cat <<USAGE
Usage:
  $0 --merged merged.vcf.gz --dna-bam DNA_TUMOR.bam --rna-bam RNA_TUMOR.bam \
     --fasta ref.fa.gz --dict ref.dict \
     --out-genotyped out.genotyped.vcf.gz --out-phased out.phased.vcf.gz [--threads 8]

Optional sample names (must match merged VCF sample columns):
  --normal-name  (default: DNA_NORMAL)
  --dna-name     (default: DNA_TUMOR)
  --rna-name     (default: RNA_TUMOR)
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
    --rna-bam) rna_bam="$2"; shift 2;;
    --fasta) fasta="$2"; shift 2;;
    --dict) dict="$2"; shift 2;;
    --out-genotyped) out_genotyped="$2"; shift 2;;
    --out-phased) out_phased="$2"; shift 2;;
    --threads) threads="$2"; shift 2;;
    --normal-name) normal_name="$2"; shift 2;;
    --dna-name) dna_name="$2"; shift 2;;
    --rna-name) rna_name="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -n "$merged_vcf" ]] || die "Missing --merged"
[[ -n "$dna_bam" ]] || die "Missing --dna-bam"
[[ -n "$rna_bam" ]] || die "Missing --rna-bam"
[[ -n "$fasta" ]] || die "Missing --fasta"
[[ -n "$dict" ]] || die "Missing --dict"
[[ -n "$out_genotyped" ]] || die "Missing --out-genotyped"
[[ -n "$out_phased" ]] || die "Missing --out-phased"

[[ -f "$merged_vcf" ]] || die "Missing merged VCF: $merged_vcf"
[[ -f "$fasta" ]] || die "Missing FASTA: $fasta"
[[ -f "${fasta}.fai" ]] || die "Missing FASTA index: ${fasta}.fai"
[[ -f "$dict" ]] || die "Missing dict: $dict"
[[ -f "$dna_bam" ]] || die "Missing DNA tumor BAM: $dna_bam"
[[ -f "$rna_bam" ]] || die "Missing RNA tumor BAM: $rna_bam"

command -v whatshap >/dev/null 2>&1 || die "whatshap not found in PATH"
command -v gatk >/dev/null 2>&1 || die "gatk not found in PATH"
command -v bcftools >/dev/null 2>&1 || die "bcftools not found in PATH"
command -v bgzip >/dev/null 2>&1 || die "bgzip not found in PATH"

bcftools index -f -t "$merged_vcf" >/dev/null 2>&1 || true

mapfile -t samples_arr < <(bcftools query -l "$merged_vcf")
samples=$(printf '%s ' "${samples_arr[@]}")
echo "[info] merged samples: $samples"

normal_in=$(resolve_sample_name "$normal_name" "${samples_arr[@]}") || die "$normal_name not found (flexible match) in merged VCF"
dna_in=$(resolve_sample_name "$dna_name" "${samples_arr[@]}") || die "$dna_name not found (flexible match) in merged VCF"
rna_in=$(resolve_sample_name "$rna_name" "${samples_arr[@]}") || die "$rna_name not found (flexible match) in merged VCF"

tmpdir=$(mktemp -d)
trap 'rm -rf "'$tmpdir'"' EXIT

# 1) Union-alleles VCF (same rows as merged, without genotypes)
union_sites_vcf="${tmpdir}/union.sites.vcf.gz"
bcftools view -G -Oz -o "$union_sites_vcf" "$merged_vcf"
bcftools index -f -t "$union_sites_vcf"

# BED of positions from merged VCF to force GATK evaluation exactly at these loci
intervals_bed="${tmpdir}/union.sites.bed"
bcftools query -f '%CHROM\t%POS\n' "$union_sites_vcf" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' > "$intervals_bed"

# 2) Genotype those same alleles in each BAM
hc_genotype_union () {
  local bam="$1"
  local out="$2"

  gatk --java-options "-Xmx16g" HaplotypeCaller -R "$fasta" -I "$bam" -L "$intervals_bed" --alleles "$union_sites_vcf" --output-mode EMIT_ALL_ACTIVE_SITES --active-probability-threshold 0.0 --native-pair-hmm-threads "$threads" -O "$out"
}

dna_raw="${tmpdir}/dna.raw.vcf.gz"
rna_raw="${tmpdir}/rna.raw.vcf.gz"

echo "[info] Genotyping union alleles in ${dna_name} BAM..."
hc_genotype_union "$dna_bam" "$dna_raw"
bcftools index -f -t "$dna_raw" >/dev/null 2>&1 || true

echo "[info] Genotyping union alleles in ${rna_name} BAM..."
hc_genotype_union "$rna_bam" "$rna_raw"
bcftools index -f -t "$rna_raw" >/dev/null 2>&1 || true

# Force sample names (in case BAM SM is different)
dna_fix="${tmpdir}/${dna_name}.vcf.gz"
rna_fix="${tmpdir}/${rna_name}.vcf.gz"

echo -e "$(bcftools query -l "$dna_raw" | head -n1)\t${dna_name}" > "${tmpdir}/dna.map"
bcftools reheader -s "${tmpdir}/dna.map" -o "$dna_fix" "$dna_raw"
bcftools index -f -t "$dna_fix"

echo -e "$(bcftools query -l "$rna_raw" | head -n1)\t${rna_name}" > "${tmpdir}/rna.map"
bcftools reheader -s "${tmpdir}/rna.map" -o "$rna_fix" "$rna_raw"
bcftools index -f -t "$rna_fix"

# Keep DNA normal genotypes from merged VCF
normal_vcf="${tmpdir}/${normal_name}.vcf.gz"
bcftools view -s "$normal_in" -Oz -o "$normal_vcf" "$merged_vcf"
bcftools index -f -t "$normal_vcf"

# Merge back to 3-sample VCF (same positions)
info_rules="KNOWN_RNAEDIT_DB:join,SOURCE_SET:join,EDIT_SIG:join"
bcftools merge --threads "$threads" -m all --info-rules "$info_rules" -Oz -o "$out_genotyped" "$normal_vcf" "$dna_fix" "$rna_fix"
bcftools index -f -t "$out_genotyped"

echo "[info] Wrote genotyped: $out_genotyped"

# 3) Read-backed phasing
pass1="${tmpdir}/pass1.${dna_name}.vcf.gz"

echo "[info] Phasing ${dna_name} (read-backed)..."
whatshap phase --reference "$fasta" --sample "$dna_name" --ignore-read-groups --output /dev/stdout "$out_genotyped" "$dna_bam" | bgzip -c > "$pass1"
bcftools index -f -t "$pass1"

echo "[info] Phasing ${rna_name} (read-backed)..."
whatshap phase --reference "$fasta" --sample "$rna_name" --ignore-read-groups --output /dev/stdout "$pass1" "$rna_bam" | bgzip -c > "$out_phased"
bcftools index -f -t "$out_phased"

echo "[info] Wrote phased: $out_phased"
