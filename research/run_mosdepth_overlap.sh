#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_mosdepth_overlap.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [--depth-threshold N]
USAGE
}

config=""
sample=""
outdir=""
depth_threshold="10"

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --depth-threshold) depth_threshold="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$config" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
if [ -z "${mosdepthdir:-}" ]; then
  if [ -n "${datadir:-}" ]; then
    mosdepthdir="${datadir%/}/reports/mosdepth"
  else
    echo "ERROR: CONFIG must define mosdepthdir (or datadir for fallback)" >&2
    exit 1
  fi
fi

if [ -z "$outdir" ]; then
  outdir="${mosdepthdir%/}/overlap_plots"
fi
mkdir -p "$outdir"

dna_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
dna_tumor="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
rna_tumor="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"

python3 "$(dirname "$0")/plot_mosdepth_overlap.py" \
  --samples "$samples" \
  --mosdepth-dir "$mosdepthdir" \
  --dna-normal-label "$dna_normal" \
  --dna-tumor-label "$dna_tumor" \
  --rna-tumor-label "$rna_tumor" \
  --depth-threshold "$depth_threshold" \
  --outdir "$outdir" \
  ${sample:+--patient "$sample"}
