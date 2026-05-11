#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_gather_mosdepth_output.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [-f] [--skip-running]
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
if [ -z "${mosdepthdir:-}" ]; then
  if [ -n "${datadir:-}" ]; then
    mosdepthdir="${datadir%/}/reports/mosdepth"
  else
    echo "ERROR: CONFIG must define mosdepthdir (or datadir for fallback)" >&2
    exit 1
  fi
fi

if [ -z "$outdir" ]; then
  outdir="${mosdepthdir%/}/gathered"
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

outfile="${outdir}/cohort.mosdepth_total_coverage.tsv"
tag="cohort"
if [ -n "$sample" ]; then
  outfile="${outdir}/${sample}.mosdepth_total_coverage.tsv"
  tag="$sample"
fi

if [ "$force" != "1" ] && [ -s "$outfile" ]; then
  echo "[skip] output already exists: ${outfile} (use -f to overwrite)"
  exit 0
fi

logroot="${outdir}/gather_mosdepth_output.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_gather_mosdepth_output"
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

python3 "${script_dir}/gather_mosdepth_output.py" \\
  --samples "${samples}" \\
  --mosdepth-dir "${mosdepthdir}" \\
  --dna-normal-label "${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}" \\
  --dna-tumor-label "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" \\
  --rna-tumor-label "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}" \\
  --outfile "${outfile}" \\
  ${sample:+--patient "${sample}"}

if [ ! -s "${outfile}" ]; then
  echo "ERROR: mosdepth coverage output missing/empty: ${outfile}" >&2
  exit 2
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
echo "[info] mosdepth coverage output file: ${outfile}"
echo ".. logs and reports saved in ${logroot}"
