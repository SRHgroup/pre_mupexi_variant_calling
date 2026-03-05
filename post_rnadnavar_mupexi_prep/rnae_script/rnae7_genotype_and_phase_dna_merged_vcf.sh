#!/usr/bin/bash
set -euo pipefail

if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  source "$PIPELINE_DEFAULTS"
fi
module load ${modules_dna_only_phase:-ngs tools htslib/1.23 bcftools/1.23 gatk/4.5.0.0 anaconda3/2025.06-1}

threads=8
hc_assembly_region_padding="${HC_ASSEMBLY_REGION_PADDING:-150}"
genotyper_mode="${dna_only_genotyper_mode:-${DNA_ONLY_GENOTYPER_MODE:-safe}}"
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
command -v bcftools >/dev/null 2>&1 || die "bcftools not found"
command -v bgzip >/dev/null 2>&1 || die "bgzip not found"
if [ "$genotyper_mode" != "safe" ]; then
  command -v gatk >/dev/null 2>&1 || die "gatk not found (required for dna_only_genotyper_mode=${genotyper_mode})"
fi

bcftools index -f -t "$merged_vcf" >/dev/null 2>&1 || true
mapfile -t samples_arr < <(bcftools query -l "$merged_vcf")
[[ ${#samples_arr[@]} -ge 2 ]] || die "Merged VCF must contain at least normal + tumor samples"

normal_in=$(resolve_sample_name "$normal_name" "${samples_arr[@]}") || die "$normal_name not found in merged VCF"
dna_in=$(resolve_sample_name "$dna_name" "${samples_arr[@]}") || die "$dna_name not found in merged VCF"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

base_genotyped="$merged_vcf"
if [ "$genotyper_mode" = "safe" ]; then
  echo "[warn] dna_only_genotyper_mode=safe: skipping HaplotypeCaller and composing TUMOR from merged VCF"
elif [ "$genotyper_mode" = "hc" ] || [ "$genotyper_mode" = "auto" ]; then
  union_sites_vcf="${tmpdir}/union.sites.vcf.gz"
  bcftools view -G -Oz -o "$union_sites_vcf" "$merged_vcf"
  bcftools index -f -t "$union_sites_vcf"

  dna_raw="${tmpdir}/dna.raw.vcf.gz"
  hc_ok=1
  if ! gatk --java-options "-Xmx16g" HaplotypeCaller \
      -R "$fasta" -I "$dna_bam" -L "$union_sites_vcf" \
      --assembly-region-padding "$hc_assembly_region_padding" \
      --alleles "$union_sites_vcf" \
      --output-mode EMIT_ALL_ACTIVE_SITES \
      --disable-tool-default-annotations \
      --active-probability-threshold 0.0 \
      --annotations-to-exclude TandemRepeat \
      --native-pair-hmm-threads "$threads" \
      -O "$dna_raw"; then
    hc_ok=0
  fi

  if [ "$hc_ok" -eq 1 ]; then
    bcftools index -f -t "$dna_raw" >/dev/null 2>&1 || true

    dna_fix="${tmpdir}/${dna_name}.vcf.gz"
    echo -e "$(bcftools query -l "$dna_raw" | head -n1)\t${dna_name}" > "${tmpdir}/dna.map"
    bcftools reheader -s "${tmpdir}/dna.map" -o "$dna_fix" "$dna_raw"
    bcftools index -f -t "$dna_fix"

    normal_vcf="${tmpdir}/${normal_name}.vcf.gz"
    bcftools view -s "$normal_in" -Oz -o "$normal_vcf" "$merged_vcf"
    bcftools index -f -t "$normal_vcf"

    info_rules="SOURCE_SET:join"
    base_genotyped="${tmpdir}/genotyped.raw.vcf.gz"
    bcftools merge --threads "$threads" -m all --info-rules "$info_rules" -Oz -o "$base_genotyped" "$normal_vcf" "$dna_fix"
    bcftools index -f -t "$base_genotyped"
  else
    if [ "$genotyper_mode" = "hc" ]; then
      die "HaplotypeCaller failed in dna_only_genotyper_mode=hc"
    fi
    echo "[warn] HaplotypeCaller failed; falling back to safe composition from merged VCF"
  fi
else
  die "Unknown dna_only_genotyper_mode: ${genotyper_mode} (allowed: auto|hc|safe)"
fi

python3 - "$base_genotyped" "$merged_vcf" "$out_genotyped" "$normal_name" "$dna_name" "$normal_in" "$dna_in" <<'PY'
import gzip
import sys

base, merged, out, normal_out, tumor_out, normal_in, tumor_in = sys.argv[1:]

def op(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def parse_header(path):
    fh = op(path, "rt")
    header = []
    samples = []
    for ln in fh:
        header.append(ln)
        if ln.startswith("#CHROM"):
            samples = ln.rstrip("\n").split("\t")[9:]
            break
    return fh, header, samples

def key(cols):
    return cols[0], int(cols[1]), cols[3], cols[4]

def cmpk(a, b):
    if a[0] != b[0]:
        return -1 if a[0] < b[0] else 1
    if a[1] != b[1]:
        return -1 if a[1] < b[1] else 1
    if a[2] != b[2]:
        return -1 if a[2] < b[2] else 1
    if a[3] != b[3]:
        return -1 if a[3] < b[3] else 1
    return 0

def next_var(fh):
    for ln in fh:
        if not ln.startswith("#"):
            return ln
    return None

def gt_missing(sample_value, fmt):
    keys = fmt.split(":")
    try:
        i = keys.index("GT")
    except ValueError:
        return True
    vals = sample_value.split(":")
    gt = vals[i] if i < len(vals) else "."
    return gt in ("", ".", "./.", ".|.")

def make_homref(fmt):
    keys = fmt.split(":")
    vals = ["."] * len(keys)
    try:
        i = keys.index("GT")
    except ValueError:
        return ":".join(vals)
    vals[i] = "0/0"
    return ":".join(vals)

bfh, bhdr, bsamples = parse_header(base)
mfh, _mhdr, msamples = parse_header(merged)

base_idx = {
    "normal": 9 + bsamples.index(normal_out) if normal_out in bsamples else (9 + bsamples.index(normal_in) if normal_in in bsamples else None),
    "tumor": 9 + bsamples.index(tumor_out) if tumor_out in bsamples else (9 + bsamples.index(tumor_in) if tumor_in in bsamples else None),
}
merged_idx = {
    "normal": 9 + msamples.index(normal_in),
    "tumor": 9 + msamples.index(tumor_in),
}
if base_idx["normal"] is None or base_idx["tumor"] is None:
    raise SystemExit("base genotyped VCF missing required samples")

mline = next_var(mfh)
filled = 0
homref_fallback = 0
records = 0
with op(out, "wt") as ofh:
    for h in bhdr:
        if h.startswith("#CHROM"):
            cols = h.rstrip("\n").split("\t")
            out_h = cols[:9] + [normal_out, tumor_out]
            ofh.write("\t".join(out_h) + "\n")
        else:
            ofh.write(h)

    while True:
        bline = next_var(bfh)
        if bline is None:
            break
        records += 1
        bc = bline.rstrip("\n").split("\t")
        bk = key(bc)

        mc = None
        while mline is not None:
            cand = mline.rstrip("\n").split("\t")
            ck = key(cand)
            c = cmpk(ck, bk)
            if c < 0:
                mline = next_var(mfh)
                continue
            if c == 0:
                mc = cand
            break

        fmt = bc[8]
        norm_b = bc[base_idx["normal"]]
        tum_b = bc[base_idx["tumor"]]
        norm_m = "./."
        tum_m = "./."
        if mc is not None:
            norm_m = mc[merged_idx["normal"]]
            tum_m = mc[merged_idx["tumor"]]

        tumor = tum_b
        if gt_missing(tumor, fmt):
            if not gt_missing(tum_m, fmt):
                tumor = tum_m
                filled += 1
            elif not gt_missing(norm_b, fmt):
                tumor = norm_b
                filled += 1
            elif not gt_missing(norm_m, fmt):
                tumor = norm_m
                filled += 1
            else:
                tumor = make_homref(fmt)
                homref_fallback += 1

        normal = norm_b
        if gt_missing(normal, fmt) and not gt_missing(norm_m, fmt):
            normal = norm_m

        out_cols = bc[:9] + [normal, tumor]
        ofh.write("\t".join(out_cols) + "\n")

print(f"[info] dna-only compose: records={records} filled_tumor={filled} homref_fallback={homref_fallback}", file=sys.stderr)
PY

bcftools +fill-tags "$out_genotyped" -Oz -o "${tmpdir}/genotyped.filled.vcf.gz" -- -t AN,AC,AF
mv "${tmpdir}/genotyped.filled.vcf.gz" "$out_genotyped"
bcftools index -f -t "$out_genotyped"

missing_tumor_gt=$(
  bcftools query -s "$dna_name" -f '[%GT\n]' "$out_genotyped" \
    | awk '($1=="./." || $1==".|." || $1=="." || $1==""){n++} END{print n+0}'
)
if [ "$missing_tumor_gt" -gt 0 ]; then
  echo "[error] dna-only validation failed: ${missing_tumor_gt} variants have missing ${dna_name} genotype in ${out_genotyped}" >&2
  exit 12
fi
echo "[info] dna-only validation passed: no missing ${dna_name} genotypes in ${out_genotyped}"

echo "[info] Wrote genotyped DNA-only VCF: $out_genotyped"

whatshap phase --reference "$fasta" --sample "$dna_name" --ignore-read-groups \
  --output /dev/stdout "$out_genotyped" "$dna_bam" | bgzip -c > "$out_phased"
bcftools index -f -t "$out_phased"

echo "[info] Wrote phased DNA-only VCF: $out_phased"
