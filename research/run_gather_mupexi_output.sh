#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_gather_mupexi_output.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [-f] [--skip-running]
USAGE
}

config=""
sample=""
outdir=""
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    -f|--force) force=1; shift ;;
    --skip-running) skip_running=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$config" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  source "$PIPELINE_DEFAULTS"
fi
source "$config"

: "${samples:?CONFIG must define samples}"
: "${mupexi_outdir:?CONFIG must define mupexi_outdir}"

if [ -z "$outdir" ]; then
  outdir="${mupexi_outdir%/}/gathered"
fi
mkdir -p "$outdir"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label:-DNA_NORMAL}"
    "${dna_tumor_label:-DNA_TUMOR}"
    "${rna_tumor_label:-RNA_TUMOR}"
    "${out_dna_normal_label:-DNA_NORMAL}"
    "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
    "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do
    value="${value%_${label}}"
  done
  printf '%s\n' "$value"
}

seen=""
snv_inputs=()
fus_inputs=()
patient_tag="cohort"
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac
  sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
  case "${sid,,}" in sample|sample_id|patient|patient_id) continue ;; esac
  patient="$(sample_base_name "$sid")"
  [ -n "$patient" ] || continue
  if printf '%s\n' "$seen" | grep -Fxq "$patient"; then
    continue
  fi
  if [ -n "$sample" ] && [ "$sample" != "$sid" ] && [ "$sample" != "$patient" ]; then
    continue
  fi
  seen="${seen}\n${patient}"
  if [ -n "$sample" ]; then
    patient_tag="$patient"
  fi
  snv="${mupexi_outdir%/}/${patient}_snv.mupexi"
  fus="${mupexi_outdir%/}/${patient}_fus.mupexi"
  if [ -f "$snv" ]; then
    snv_inputs+=("--snv-input" "${patient}=${snv}")
  fi
  if [ -f "$fus" ]; then
    fus_inputs+=("--fus-input" "${patient}=${fus}")
  fi
done < "$samples"

snv_out="${outdir}/${patient_tag}.snv.mupexi.tsv"
fus_out="${outdir}/${patient_tag}.fus.mupexi.tsv"
if [ "$force" != "1" ] && [ -s "$snv_out" ] && [ -s "$fus_out" ]; then
  echo "[skip] outputs already exist: $snv_out and $fus_out (use -f to overwrite)"
  exit 0
fi

python3 "$(cd "$(dirname "$0")/.." && pwd)/research/gather_mupexi_output.py" \
  "${snv_inputs[@]}" \
  "${fus_inputs[@]}" \
  --snv-outfile "$snv_out" \
  --fus-outfile "$fus_out"

echo "[done] SNV gathered -> $snv_out"
echo "[done] FUS gathered -> $fus_out"
