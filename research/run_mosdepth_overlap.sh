#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_mosdepth_overlap.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [--depth-threshold N] [-f] [--skip-running]
USAGE
}

config=""
sample=""
outdir=""
depth_threshold="10"
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --depth-threshold) depth_threshold="$2"; shift 2 ;;
    -f|--force) force=1; shift ;;
    --skip-running) skip_running=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$config" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

# Optional shared toolchain defaults (set by run_pipeline wrappers)
if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "$PIPELINE_DEFAULTS"
fi

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
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

dna_normal="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
dna_tumor="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
rna_tumor="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

summary_tsv="${outdir}/mosdepth_overlap_summary.tsv"
categories_tsv="${outdir}/mosdepth_overlap_categories.tsv"
missing_tsv="${outdir}/mosdepth_overlap_missing_inputs.tsv"
stacked_png="${outdir}/mosdepth_overlap_stacked.png"
stacked_svg="${outdir}/mosdepth_overlap_stacked.svg"
jaccard_png="${outdir}/mosdepth_overlap_jaccard_heatmap.png"
jaccard_svg="${outdir}/mosdepth_overlap_jaccard_heatmap.svg"

if [ "$force" != "1" ] && [ -s "$summary_tsv" ] && [ -s "$categories_tsv" ] && [ -s "$stacked_png" ] && [ -s "$stacked_svg" ] && [ -s "$jaccard_png" ] && [ -s "$jaccard_svg" ]; then
  echo "[skip] mosdepth overlap outputs already exist in ${outdir} (use -f to overwrite)"
  exit 0
fi

logroot="${outdir}/mosdepth_overlap.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_mosdepth_overlap"
tag="cohort"
[ -n "$sample" ] && tag="$sample"
marker="${logdir}/submitted.${prefix}.${tag}.jobid"

pbs_state_for_jobid() {
  local jid="$1"
  local line
  line="$(qstat -f "$jid" 2>/dev/null | awk -F' = ' '/job_state =/{print $2; exit}' || true)"
  case "$line" in
    R|E) printf '%s\n' "RUNNING" ;;
    Q|H|W|T|S) printf '%s\n' "QUEUED" ;;
    *) printf '%s\n' "" ;;
  esac
}

if [ "$skip_running" = "1" ] && [ -f "$marker" ]; then
  prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
  if [ -n "$prev_jobid" ]; then
    st="$(pbs_state_for_jobid "$prev_jobid")"
    if [ "$st" = "RUNNING" ] || [ "$st" = "QUEUED" ]; then
      echo "[skip-running] ${prefix}.${tag}: active job ${prev_jobid} (${st})"
      exit 0
    fi
  fi
fi

active_jobid=""
if command -v qselect >/dev/null 2>&1; then
  active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "${prefix}.${tag}" 2>/dev/null | head -n1 || true)"
fi
if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
  active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="${prefix}.${tag}" '$4==n {print $1; exit}')"
fi
if [ -n "$active_jobid" ]; then
  echo "[skip] ${prefix}.${tag}: scheduler already has active job ${active_jobid}"
  exit 0
fi

runscript="${logdir}/run.${tag}.${prefix}.sh"
cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${script_dir}/plot_mosdepth_overlap.py" \\
  --samples "${samples}" \\
  --mosdepth-dir "${mosdepthdir}" \\
  --dna-normal-label "${dna_normal}" \\
  --dna-tumor-label "${dna_tumor}" \\
  --rna-tumor-label "${rna_tumor}" \\
  --depth-threshold "${depth_threshold}" \\
  --outdir "${outdir}" \\
  ${sample:+--patient "${sample}"}

if [ ! -s "${summary_tsv}" ]; then
  echo "ERROR: summary output missing/empty: ${summary_tsv}" >&2
  exit 2
fi
if [ ! -s "${categories_tsv}" ]; then
  echo "ERROR: category output missing/empty: ${categories_tsv}" >&2
  exit 3
fi
SCRIPT
chmod +x "$runscript"

qsub_opts=()
[ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
[ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
qsub_opts+=(-N "${prefix}.${tag}")
qsub_opts+=(-o "${repdir}/${prefix}.${tag}.o\$PBS_JOBID")
qsub_opts+=(-e "${repdir}/${prefix}.${tag}.e\$PBS_JOBID")

jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
printf '%s\n' "$jobid" > "$marker"
echo "[submit] ${prefix}.${tag}: jobid=${jobid}"
echo "[info] outputs target dir: ${outdir}"
echo ".. logs and reports saved in ${logroot}"
