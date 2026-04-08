#!/usr/bin/bash
set -euo pipefail

module load ngs tools htslib/1.23 bcftools/1.23 anaconda3/2025.06-1

threads=4
rna_keep_all=0
rna_known_regex='REDIportal|RADAR|Asaoka|APOBEC3_MOTIF'   # match within INFO/KNOWN_RNAEDIT_DB

usage() {
  cat <<EOF
Usage:
  $0 -g GERMLINE.vcf.gz -d DNA_TUMOR.vcf.gz -r RNA_TUMOR.vcf.gz -o merged.vcf.gz [options]

Required:
  -g   germline VCF (bgzipped; HaplotypeCaller-style OK)
  -d   DNA tumor somatic VCF (bgzipped; Mutect2-style)
  -r   RNA tumor VCF (bgzipped; Mutect2-style; should have KNOWN_RNAEDIT_DB if using known-only)
  -o   output merged VCF (MUST end with .vcf.gz)

Options:
  -t INT                 threads (default: 4)
  --rna-keep-all         do NOT filter RNA to known sites (default: filter to known)
  --rna-known-regex STR  regex used on INFO/KNOWN_RNAEDIT_DB (default: 'REDIportal|RADAR|Asaoka|APOBEC3_MOTIF')

Notes:
  - Output sample names are EXACTLY: DNA_NORMAL, DNA_TUMOR, RNA_TUMOR
  - Adds INFO flags: GERMLINE on germline records; SOMATIC on DNA tumor records
  - Strips allele-specific INFO tags (AS_*) from Mutect2 VCFs to avoid merge crashes
  - Preserves RNA annotations at overlapping sites by using mergeable INFO strings:
      SOURCE_SET (GERMLINE/SOMATIC/RNA_EDIT) and EDIT_SIG (ADAR/APOBEC3)
    and restores flag tags after merge.
EOF
}

germ_vcf=""
dna_vcf=""
rna_vcf=""
out_vcf=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) germ_vcf="$2"; shift 2;;
    -d) dna_vcf="$2"; shift 2;;
    -r) rna_vcf="$2"; shift 2;;
    -o) out_vcf="$2"; shift 2;;
    -t) threads="$2"; shift 2;;
    --rna-keep-all) rna_keep_all=1; shift 1;;
    --rna-known-regex) rna_known_regex="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

if [[ -z "${germ_vcf}" || -z "${dna_vcf}" || -z "${rna_vcf}" || -z "${out_vcf}" ]]; then
  usage >&2
  exit 1
fi

if [[ "${out_vcf}" != *.vcf.gz ]]; then
  echo "ERROR: -o must end with .vcf.gz (bgzipped). You gave: ${out_vcf}" >&2
  exit 1
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

count_records() { bcftools view -H "$1" | wc -l; }
count_flag() { bcftools view -H -i "INFO/$1" "$2" | wc -l; }
count_str_notempty() { bcftools view -H -i "INFO/$1!=\"\" && INFO/$1!=\".\"" "$2" | wc -l; }

# -------------------------
# helper: strip AS_* INFO fields (Mutect2 allele-specific tags)
# -------------------------
strip_as_info() {
  local in_vcf="$1"
  local out_vcf="$2"

  local as_tags
  as_tags=$(bcftools view -h "${in_vcf}" \
    | sed -n 's/^##INFO=<ID=\(AS_[^,>]*\).*/INFO\/\1/p' \
    | paste -sd, -)

  if [[ -n "${as_tags}" ]]; then
    echo "[info] Stripping allele-specific INFO tags (AS_*) from: ${in_vcf}"
    echo "[info] Removing: ${as_tags}"
    bcftools annotate --threads "${threads}" -x "${as_tags}" -Oz -o "${out_vcf}" "${in_vcf}"
    bcftools index -t "${out_vcf}"
  else
    cp -f "${in_vcf}" "${out_vcf}"
    if [[ -f "${in_vcf}.tbi" ]]; then
      cp -f "${in_vcf}.tbi" "${out_vcf}.tbi"
    else
      bcftools index -t "${out_vcf}" >/dev/null 2>&1 || true
    fi
  fi
}

# -------------------------
# helper: add INFO flag to all records (GERMLINE/SOMATIC)
# -------------------------
add_info_flag_py="${tmpdir}/add_info_flag.py"
cat > "${add_info_flag_py}" <<'PY'
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
chmod +x "${add_info_flag_py}"

# -------------------------
# helper: add/append INFO string key=value (SOURCE_SET, etc.)
# (stores comma-separated unique values)
# -------------------------
add_info_list_py="${tmpdir}/add_info_list.py"
cat > "${add_info_list_py}" <<'PY'
#!/usr/bin/env python3
import sys

key = sys.argv[1]
value = sys.argv[2]
desc = sys.argv[3]

hdr_line = f'##INFO=<ID={key},Number=.,Type=String,Description="{desc}">\n'

seen_hdr = False
for line in sys.stdin:
    if line.startswith("##INFO=<ID=" + key + ","):
        seen_hdr = True
        sys.stdout.write(line)
        continue

    if line.startswith("##"):
        sys.stdout.write(line)
        continue

    if line.startswith("#CHROM"):
        if not seen_hdr:
            sys.stdout.write(hdr_line)
        sys.stdout.write(line)
        continue

    if not line or line[0] == "#":
        sys.stdout.write(line)
        continue

    toks = line.rstrip("\n").split("\t")
    if len(toks) < 8:
        continue

    info = toks[7]
    fields = [] if info == "." else info.split(";")

    found = False
    for i,f in enumerate(fields):
        if f.startswith(key + "="):
            found = True
            cur = f.split("=",1)[1]
            cur_vals = [x for x in cur.split(",") if x]
            if value not in cur_vals:
                cur_vals.append(value)
            fields[i] = key + "=" + ",".join(cur_vals)
            break
    if not found:
        fields.append(key + "=" + value)

    toks[7] = ";".join([x for x in fields if x]) if fields else "."
    sys.stdout.write("\t".join(toks) + "\n")
PY
chmod +x "${add_info_list_py}"

# -------------------------
# helper: in RNA VCF, convert ADAR_SIG/APOBEC3_SIG flags into mergeable string INFO EDIT_SIG
# -------------------------
rna_make_edit_sig_py="${tmpdir}/rna_make_edit_sig.py"
cat > "${rna_make_edit_sig_py}" <<'PY'
#!/usr/bin/env python3
import sys

key = "EDIT_SIG"
desc = 'RNA editing signature label(s) derived from flags ADAR_SIG/APOBEC3_SIG (values: ADAR,APOBEC3)'

hdr_line = f'##INFO=<ID={key},Number=.,Type=String,Description="{desc}">\n'
seen_hdr = False

for line in sys.stdin:
    if line.startswith("##INFO=<ID=" + key + ","):
        seen_hdr = True
        sys.stdout.write(line)
        continue
    if line.startswith("##"):
        sys.stdout.write(line)
        continue
    if line.startswith("#CHROM"):
        if not seen_hdr:
            sys.stdout.write(hdr_line)
        sys.stdout.write(line)
        continue
    if not line or line[0] == "#":
        sys.stdout.write(line)
        continue

    toks = line.rstrip("\n").split("\t")
    if len(toks) < 8:
        continue

    info = toks[7]
    fields = [] if info == "." else info.split(";")
    field_set = set(fields)

    sigs = []
    if "ADAR_SIG" in field_set:
        sigs.append("ADAR")
    if "APOBEC3_SIG" in field_set:
        sigs.append("APOBEC3")

    # If no signatures, leave as-is
    if not sigs:
        sys.stdout.write(line)
        continue

    # Add/append EDIT_SIG
    found = False
    for i,f in enumerate(fields):
        if f.startswith(key + "="):
            found = True
            cur = f.split("=",1)[1]
            cur_vals = [x for x in cur.split(",") if x]
            for s in sigs:
                if s not in cur_vals:
                    cur_vals.append(s)
            fields[i] = key + "=" + ",".join(cur_vals)
            break
    if not found:
        fields.append(key + "=" + ",".join(sigs))

    toks[7] = ";".join([x for x in fields if x]) if fields else "."
    sys.stdout.write("\t".join(toks) + "\n")
PY
chmod +x "${rna_make_edit_sig_py}"

# -------------------------
# helper: after merge, restore flags from SOURCE_SET/EDIT_SIG
# -------------------------
restore_flags_py="${tmpdir}/restore_flags.py"
cat > "${restore_flags_py}" <<'PY'
#!/usr/bin/env python3
import sys

# Ensure headers exist (won't duplicate if already present)
hdr_needed = [
  ('GERMLINE','Variant present in germline callset'),
  ('SOMATIC','Variant present in DNA tumor somatic callset'),
  ('ADAR_SIG','RNA editing signature: ADAR'),
  ('APOBEC3_SIG','RNA editing signature: APOBEC3'),
]
seen = {k: False for k,_ in hdr_needed}

def emit_missing_headers():
    for k,desc in hdr_needed:
        if not seen[k]:
            sys.stdout.write(f'##INFO=<ID={k},Number=0,Type=Flag,Description="{desc}">\n')

for line in sys.stdin:
    if line.startswith("##INFO=<ID="):
        for k,_ in hdr_needed:
            if line.startswith(f"##INFO=<ID={k},"):
                seen[k] = True
        sys.stdout.write(line)
        continue

    if line.startswith("##"):
        sys.stdout.write(line)
        continue

    if line.startswith("#CHROM"):
        emit_missing_headers()
        sys.stdout.write(line)
        continue

    if not line or line[0] == "#":
        sys.stdout.write(line)
        continue

    toks = line.rstrip("\n").split("\t")
    if len(toks) < 8:
        continue

    info = toks[7]
    fields = [] if info == "." else info.split(";")
    fset = set(fields)

    # parse SOURCE_SET and EDIT_SIG values (if present)
    source_vals = set()
    edit_vals = set()
    for f in fields:
        if f.startswith("SOURCE_SET="):
            source_vals |= set([x for x in f.split("=",1)[1].split(",") if x])
        if f.startswith("EDIT_SIG="):
            edit_vals |= set([x for x in f.split("=",1)[1].split(",") if x])

    # restore flags
    if "GERMLINE" in source_vals and "GERMLINE" not in fset:
        fields.append("GERMLINE")
    if "SOMATIC" in source_vals and "SOMATIC" not in fset:
        fields.append("SOMATIC")
    if "ADAR" in edit_vals and "ADAR_SIG" not in fset:
        fields.append("ADAR_SIG")
    if "APOBEC3" in edit_vals and "APOBEC3_SIG" not in fset:
        fields.append("APOBEC3_SIG")

    toks[7] = ";".join([x for x in fields if x]) if fields else "."
    sys.stdout.write("\t".join(toks) + "\n")
PY
chmod +x "${restore_flags_py}"

# -------------------------
# sample picking
# -------------------------
pick_first_sample() { bcftools query -l "$1" | head -n 1; }

pick_tumorish_sample() {
  local vcf="$1"
  local samples
  mapfile -t samples < <(bcftools query -l "${vcf}")
  [[ "${#samples[@]}" -eq 0 ]] && { echo ""; return; }
  for s in "${samples[@]}"; do
    if [[ "${s^^}" =~ TUMOU?R ]]; then echo "${s}"; return; fi
  done
  for s in "${samples[@]}"; do
    if [[ "${s^^}" =~ RNA ]]; then echo "${s}"; return; fi
  done
  echo "${samples[-1]}"
}

# -------------------------
# 1) Germline -> DNA_NORMAL + GERMLINE + SOURCE_SET=GERMLINE
# -------------------------
germ_sample="$(pick_first_sample "${germ_vcf}")"
[[ -z "${germ_sample}" ]] && { echo "No sample found in germline VCF: ${germ_vcf}" >&2; exit 1; }

germ_sub="${tmpdir}/germ.sub.vcf.gz"
bcftools view --threads "${threads}" -s "${germ_sample}" -Oz -o "${germ_sub}" "${germ_vcf}"
bcftools index -t "${germ_sub}"

germ_map="${tmpdir}/germ.samples.tsv"
printf "%s\tDNA_NORMAL\n" "${germ_sample}" > "${germ_map}"

germ_norm="${tmpdir}/germ.DNA_NORMAL.vcf.gz"
bcftools reheader -s "${germ_map}" -o "${germ_norm}" "${germ_sub}"
bcftools index -t "${germ_norm}"

germ_tag="${tmpdir}/germ.DNA_NORMAL.tagged.vcf.gz"
bcftools view --threads "${threads}" "${germ_norm}" \
  | python3 "${add_info_flag_py}" GERMLINE "Variant present in germline callset" \
  | python3 "${add_info_list_py}" SOURCE_SET GERMLINE "Which callset(s) this record originated from (GERMLINE,SOMATIC,RNA_EDIT)" \
  | bgzip -c > "${germ_tag}"
bcftools index -t "${germ_tag}"

# -------------------------
# 2) DNA somatic -> DNA_TUMOR + SOMATIC + SOURCE_SET=SOMATIC (strip AS_*)
# -------------------------
dna_sample="$(pick_tumorish_sample "${dna_vcf}")"
[[ -z "${dna_sample}" ]] && { echo "No tumor-ish sample found in DNA VCF: ${dna_vcf}" >&2; exit 1; }

dna_sub="${tmpdir}/dna.sub.vcf.gz"
bcftools view --threads "${threads}" -s "${dna_sample}" -Oz -o "${dna_sub}" "${dna_vcf}"
bcftools index -t "${dna_sub}"

dna_map="${tmpdir}/dna.samples.tsv"
printf "%s\tDNA_TUMOR\n" "${dna_sample}" > "${dna_map}"

dna_norm="${tmpdir}/dna.DNA_TUMOR.vcf.gz"
bcftools reheader -s "${dna_map}" -o "${dna_norm}" "${dna_sub}"
bcftools index -t "${dna_norm}"

dna_clean="${tmpdir}/dna.DNA_TUMOR.noAS.vcf.gz"
strip_as_info "${dna_norm}" "${dna_clean}"

dna_tag="${tmpdir}/dna.DNA_TUMOR.tagged.vcf.gz"
bcftools view --threads "${threads}" "${dna_clean}" \
  | python3 "${add_info_flag_py}" SOMATIC "Variant present in DNA tumor somatic callset" \
  | python3 "${add_info_list_py}" SOURCE_SET SOMATIC "Which callset(s) this record originated from (GERMLINE,SOMATIC,RNA_EDIT)" \
  | bgzip -c > "${dna_tag}"
bcftools index -t "${dna_tag}"

# -------------------------
# 3) RNA editing -> RNA_TUMOR + SOURCE_SET=RNA_EDIT + EDIT_SIG derived from flags (strip AS_*), optionally filter to known sites
# -------------------------
rna_sample="$(pick_tumorish_sample "${rna_vcf}")"
[[ -z "${rna_sample}" ]] && { echo "No tumor-ish sample found in RNA VCF: ${rna_vcf}" >&2; exit 1; }

rna_sub="${tmpdir}/rna.sub.vcf.gz"
bcftools view --threads "${threads}" -s "${rna_sample}" -Oz -o "${rna_sub}" "${rna_vcf}"
bcftools index -t "${rna_sub}"

rna_map="${tmpdir}/rna.samples.tsv"
printf "%s\tRNA_TUMOR\n" "${rna_sample}" > "${rna_map}"

rna_norm="${tmpdir}/rna.RNA_TUMOR.vcf.gz"
bcftools reheader -s "${rna_map}" -o "${rna_norm}" "${rna_sub}"
bcftools index -t "${rna_norm}"

rna_clean="${tmpdir}/rna.RNA_TUMOR.noAS.vcf.gz"
strip_as_info "${rna_norm}" "${rna_clean}"

# add SOURCE_SET=RNA_EDIT and EDIT_SIG=...
rna_annot="${tmpdir}/rna.RNA_TUMOR.annot.vcf.gz"
bcftools view --threads "${threads}" "${rna_clean}" \
  | python3 "${add_info_list_py}" SOURCE_SET RNA_EDIT "Which callset(s) this record originated from (GERMLINE,SOMATIC,RNA_EDIT)" \
  | python3 "${rna_make_edit_sig_py}" \
  | bgzip -c > "${rna_annot}"
bcftools index -t "${rna_annot}"

rna_final="${tmpdir}/rna.final.vcf.gz"
if [[ "${rna_keep_all}" -eq 1 ]]; then
  cp -f "${rna_annot}" "${rna_final}"
  cp -f "${rna_annot}.tbi" "${rna_final}.tbi"
else
  if ! bcftools view -h "${rna_annot}" | grep -q '##INFO=<ID=KNOWN_RNAEDIT_DB'; then
    echo "ERROR: RNA VCF lacks INFO/KNOWN_RNAEDIT_DB but known-only is enabled." >&2
    echo "Annotate RNA first, or rerun with --rna-keep-all." >&2
    exit 1
  fi

  bcftools view --threads "${threads}" \
    -i "INFO/KNOWN_RNAEDIT_DB!=\"\" && INFO/KNOWN_RNAEDIT_DB!=\".\" && INFO/KNOWN_RNAEDIT_DB~\"${rna_known_regex}\"" \
    -Oz -o "${rna_final}" "${rna_annot}"
  bcftools index -t "${rna_final}"

  # hard fail if you filtered everything (this is exactly what happened to you)
  if [[ "$(count_records "${rna_final}")" -eq 0 ]]; then
    echo "ERROR: RNA known-sites filter produced 0 records." >&2
    echo "This means KNOWN_RNAEDIT_DB is empty/absent for all RNA records in the file you provided." >&2
    echo "Quick check:" >&2
    echo "  bcftools query -f '%INFO/KNOWN_RNAEDIT_DB\\n' '${rna_annot}' | head" >&2
    echo "  bcftools view -H -i 'INFO/KNOWN_RNAEDIT_DB!=\"\" && INFO/KNOWN_RNAEDIT_DB!=\".\"' '${rna_annot}' | head" >&2
    exit 1
  fi
fi

# -------------------------
# Pre-merge sanity
# -------------------------
echo "[info] Pre-merge counts:"
echo "  germ_tag records: $(count_records "${germ_tag}")"
echo "  dna_tag  records: $(count_records "${dna_tag}")"
echo "  rna_final records: $(count_records "${rna_final}")"
echo "  rna_final ADAR_SIG:    $(count_flag ADAR_SIG "${rna_final}")"
echo "  rna_final APOBEC3_SIG: $(count_flag APOBEC3_SIG "${rna_final}")"
echo "  rna_final EDIT_SIG:    $(count_str_notempty EDIT_SIG "${rna_final}")"
echo "  rna_final knownsites:  $(count_str_notempty KNOWN_RNAEDIT_DB "${rna_final}")"

# -------------------------
# 4) Merge (ONLY string tags in info-rules; flags are not supported there)
# -------------------------
info_rules="KNOWN_RNAEDIT_DB:join,SOURCE_SET:join,EDIT_SIG:join"

merged_tmp="${tmpdir}/merged.tmp.vcf.gz"
bcftools merge --threads "${threads}" -m all \
  --info-rules "${info_rules}" \
  -Oz -o "${merged_tmp}" \
  "${germ_tag}" "${dna_tag}" "${rna_final}"
bcftools index -t "${merged_tmp}"

# Restore flags from SOURCE_SET/EDIT_SIG into final output
bcftools view --threads "${threads}" "${merged_tmp}" \
  | python3 "${restore_flags_py}" \
  | bgzip -c > "${out_vcf}"
bcftools index -t "${out_vcf}"

echo "Wrote: ${out_vcf}"
echo "Samples in merged:"
bcftools query -l "${out_vcf}"

echo "[info] Post-merge tag presence (record counts):"
echo "  merged GERMLINE:    $(count_flag GERMLINE "${out_vcf}")"
echo "  merged SOMATIC:     $(count_flag SOMATIC "${out_vcf}")"
echo "  merged ADAR_SIG:    $(count_flag ADAR_SIG "${out_vcf}")"
echo "  merged APOBEC3_SIG: $(count_flag APOBEC3_SIG "${out_vcf}")"
echo "  merged EDIT_SIG:    $(count_str_notempty EDIT_SIG "${out_vcf}")"
echo "  merged knownsites:  $(count_str_notempty KNOWN_RNAEDIT_DB "${out_vcf}")"
