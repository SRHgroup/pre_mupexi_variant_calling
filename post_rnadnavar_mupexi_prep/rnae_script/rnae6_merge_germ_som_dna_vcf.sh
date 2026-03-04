#!/usr/bin/bash
set -euo pipefail

module load ngs tools htslib/1.23 bcftools/1.23 anaconda3/2025.06-1

threads=4
out_normal="DNA_NORMAL"
out_dna="DNA_TUMOR"

germ_vcf=""
dna_vcf=""
out_vcf=""

usage() {
  cat <<USAGE
Usage:
  $0 -g GERMLINE.vcf.gz -d DNA_SOMATIC.vcf.gz -o merged.vcf.gz [options]

Required:
  -g   germline VCF (bgzipped)
  -d   DNA tumor-vs-normal somatic VCF (bgzipped)
  -o   output merged VCF (.vcf.gz)

Options:
  -t INT             threads (default: 4)
  --out-normal STR   output normal sample name (default: DNA_NORMAL)
  --out-dna STR      output tumor sample name (default: DNA_TUMOR)
USAGE
}

die(){ echo "[error] $*" >&2; exit 1; }

normalize_to_bgzip() {
  local in_vcf="$1"
  local out_vcfgz="$2"
  bcftools view -Oz -o "$out_vcfgz" "$in_vcf"
  bcftools index -f -t "$out_vcfgz"
}

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

strip_as_info() {
  local in_vcf="$1"
  local out_vcf="$2"

  local as_tags
  as_tags=$(bcftools view -h "${in_vcf}" | sed -n 's/^##INFO=<ID=\(AS_[^,>]*\).*/INFO\/\1/p' | paste -sd, -)

  if [[ -n "${as_tags}" ]]; then
    bcftools annotate --threads "${threads}" -x "${as_tags}" -Oz -o "${out_vcf}" "${in_vcf}"
    bcftools index -f -t "${out_vcf}"
  else
    normalize_to_bgzip "${in_vcf}" "${out_vcf}"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) germ_vcf="$2"; shift 2;;
    -d) dna_vcf="$2"; shift 2;;
    -o) out_vcf="$2"; shift 2;;
    -t) threads="$2"; shift 2;;
    --out-normal) out_normal="$2"; shift 2;;
    --out-dna) out_dna="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -n "$germ_vcf" ]] || die "Missing -g"
[[ -n "$dna_vcf" ]] || die "Missing -d"
[[ -n "$out_vcf" ]] || die "Missing -o"
[[ "$out_vcf" == *.vcf.gz ]] || die "-o must end with .vcf.gz"
[[ -f "$germ_vcf" ]] || die "Missing germline VCF: $germ_vcf"
[[ -f "$dna_vcf" ]] || die "Missing DNA somatic VCF: $dna_vcf"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

germ_bgz="${tmpdir}/germ.input.vcf.gz"
dna_bgz="${tmpdir}/dna.input.vcf.gz"
normalize_to_bgzip "$germ_vcf" "$germ_bgz"
normalize_to_bgzip "$dna_vcf" "$dna_bgz"

add_info_flag_py="${tmpdir}/add_info_flag.py"
cat > "$add_info_flag_py" <<'PY'
#!/usr/bin/env python3
import sys
flag = sys.argv[1]
desc = sys.argv[2]
seen = False
for line in sys.stdin:
    if line.startswith("##INFO=<ID=" + flag + ","):
        seen = True
        sys.stdout.write(line)
        continue
    if line.startswith("##"):
        sys.stdout.write(line)
        continue
    if line.startswith("#CHROM"):
        if not seen:
            sys.stdout.write(f'##INFO=<ID={flag},Number=0,Type=Flag,Description="{desc}">\n')
        sys.stdout.write(line)
        continue
    if not line or line[0] == "#":
        sys.stdout.write(line)
        continue
    toks = line.rstrip("\n").split("\t")
    if len(toks) < 8:
        continue
    info = toks[7]
    parts = [] if info == "." else info.split(";")
    if flag not in parts:
        parts.append(flag)
    toks[7] = ";".join([p for p in parts if p]) if parts else "."
    sys.stdout.write("\t".join(toks) + "\n")
PY
chmod +x "$add_info_flag_py"

add_info_list_py="${tmpdir}/add_info_list.py"
cat > "$add_info_list_py" <<'PY'
#!/usr/bin/env python3
import sys
key = sys.argv[1]
value = sys.argv[2]
desc = sys.argv[3]
hdr = f'##INFO=<ID={key},Number=.,Type=String,Description="{desc}">\n'
seen = False
for line in sys.stdin:
    if line.startswith("##INFO=<ID=" + key + ","):
        seen = True
        sys.stdout.write(line)
        continue
    if line.startswith("##"):
        sys.stdout.write(line)
        continue
    if line.startswith("#CHROM"):
        if not seen:
            sys.stdout.write(hdr)
        sys.stdout.write(line)
        continue
    if not line or line[0] == '#':
        sys.stdout.write(line)
        continue
    toks = line.rstrip("\n").split("\t")
    if len(toks) < 8:
        continue
    info = toks[7]
    fields = [] if info == "." else info.split(";")
    found = False
    for i, f in enumerate(fields):
        if f.startswith(key + "="):
            found = True
            cur = [x for x in f.split("=", 1)[1].split(",") if x]
            if value not in cur:
                cur.append(value)
            fields[i] = key + "=" + ",".join(cur)
            break
    if not found:
        fields.append(key + "=" + value)
    toks[7] = ";".join([x for x in fields if x]) if fields else "."
    sys.stdout.write("\t".join(toks) + "\n")
PY
chmod +x "$add_info_list_py"

mapfile -t gsamples < <(bcftools query -l "$germ_bgz")
[[ ${#gsamples[@]} -ge 1 ]] || die "No samples found in germline VCF"
germ_in="${gsamples[0]}"

mapfile -t dsamples < <(bcftools query -l "$dna_bgz")
[[ ${#dsamples[@]} -ge 1 ]] || die "No samples found in DNA somatic VCF"
dna_in="$(resolve_sample_name "$out_dna" "${dsamples[@]}")" || {
  # fallback: pick first non-normal-like sample
  dna_in=""
  for s in "${dsamples[@]}"; do
    if [[ ! "${s^^}" =~ DNA_NORMAL$ ]] && [[ ! "${s^^}" =~ _N$ ]]; then
      dna_in="$s"
      break
    fi
  done
  [[ -n "$dna_in" ]] || dna_in="${dsamples[-1]}"
}

germ_named="${tmpdir}/germ.named.vcf.gz"
echo -e "${germ_in}\t${out_normal}" > "${tmpdir}/germ.map"
germ_rehead="${tmpdir}/germ.reheader.vcf"
bcftools reheader -s "${tmpdir}/germ.map" -o "$germ_rehead" "$germ_bgz"
bcftools view -Oz -o "$germ_named" "$germ_rehead"
bcftools index -f -t "$germ_named"

germ_proc="${tmpdir}/germ.proc.vcf.gz"
bcftools view -Ou "$germ_named" \
  | python3 "$add_info_flag_py" GERMLINE "Germline call" \
  | python3 "$add_info_list_py" SOURCE_SET GERMLINE "Which callset(s) this record originated from (GERMLINE,SOMATIC)" \
  | bcftools view -Oz -o "$germ_proc"
bcftools index -f -t "$germ_proc"

dna_noas="${tmpdir}/dna.noas.vcf.gz"
strip_as_info "$dna_bgz" "$dna_noas"

dna_one="${tmpdir}/dna.one.vcf.gz"
bcftools view -s "$dna_in" -Oz -o "$dna_one" "$dna_noas"
bcftools index -f -t "$dna_one"

echo -e "${dna_in}\t${out_dna}" > "${tmpdir}/dna.map"
dna_named="${tmpdir}/dna.named.vcf.gz"
dna_rehead="${tmpdir}/dna.reheader.vcf"
bcftools reheader -s "${tmpdir}/dna.map" -o "$dna_rehead" "$dna_one"
bcftools view -Oz -o "$dna_named" "$dna_rehead"
bcftools index -f -t "$dna_named"

dna_proc="${tmpdir}/dna.proc.vcf.gz"
bcftools view -Ou "$dna_named" \
  | python3 "$add_info_flag_py" SOMATIC "Somatic DNA call" \
  | python3 "$add_info_list_py" SOURCE_SET SOMATIC "Which callset(s) this record originated from (GERMLINE,SOMATIC)" \
  | bcftools view -Oz -o "$dna_proc"
bcftools index -f -t "$dna_proc"

mkdir -p "$(dirname "$out_vcf")"
info_rules="SOURCE_SET:join"
bcftools merge --threads "$threads" -m all --info-rules "$info_rules" -Oz -o "$out_vcf" "$germ_proc" "$dna_proc"
bcftools index -f -t "$out_vcf"

echo "[info] Wrote merged DNA-only VCF: $out_vcf"
echo "[info] Records: $(bcftools view -H "$out_vcf" | wc -l)"
