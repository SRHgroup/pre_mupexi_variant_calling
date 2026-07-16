#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_gather_mupexi_by_suffix.sh -c CONFIG --suffix _neospl.mupexi [-s PATIENT] [-o OUTDIR] [--input-dir DIR] [--outfile TSV] [-f] [--skip-running]

Examples:
  bash research/run_gather_mupexi_by_suffix.sh -c CONFIG --suffix _snv.mupexi
  bash research/run_gather_mupexi_by_suffix.sh -c CONFIG --suffix _fus.mupexi
  bash research/run_gather_mupexi_by_suffix.sh -c CONFIG --suffix _neospl.mupexi
  bash research/run_gather_mupexi_by_suffix.sh -c CONFIG --suffix _neojunctions.mupexi
USAGE
}

config=""
sample=""
outdir=""
input_dir=""
outfile=""
suffix=""
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="${2:-}"; shift 2 ;;
    -s|--sample) sample="${2:-}"; shift 2 ;;
    -o|--outdir) outdir="${2:-}"; shift 2 ;;
    --input-dir) input_dir="${2:-}"; shift 2 ;;
    --outfile) outfile="${2:-}"; shift 2 ;;
    --suffix) suffix="${2:-}"; shift 2 ;;
    -f|--force) force=1; shift ;;
    --skip-running) skip_running=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[ -n "$config" ] || { usage >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }
[ -n "$suffix" ] || { echo "ERROR: --suffix is required, e.g. --suffix _neospl.mupexi" >&2; exit 1; }
case "$suffix" in
  *.mupexi) ;;
  *) echo "ERROR: --suffix must end with .mupexi, got: $suffix" >&2; exit 1 ;;
esac

if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "$PIPELINE_DEFAULTS"
fi
# shellcheck disable=SC1090
source "$config"

: "${mupexi_outdir:?CONFIG must define mupexi_outdir}"

if [ -z "$input_dir" ]; then
  input_dir="$mupexi_outdir"
fi
if [ -z "$outdir" ]; then
  outdir="${mupexi_outdir%/}/gathered"
fi

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

patient_tag="cohort"
patient_arg=()
if [ -n "$sample" ]; then
  patient_tag="$(sample_base_name "$sample")"
  patient_arg=(--patient "$patient_tag")
fi

suffix_label="${suffix#_}"
suffix_label="${suffix_label//\//_}"
suffix_label="${suffix_label// /_}"
if [ -z "$outfile" ]; then
  outfile="${outdir%/}/${patient_tag}.${suffix_label}.tsv"
fi

if [ "$force" != "1" ] && [ -s "$outfile" ]; then
  echo "[skip] output already exists: $outfile (use -f to overwrite)"
  exit 0
fi

mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

safe_suffix="${suffix//[^A-Za-z0-9_.-]/_}"
safe_suffix="${safe_suffix#_}"
prefix="research_gather_mupexi_by_suffix"
job_tag="${patient_tag}.${safe_suffix}"
logroot="${outdir}/gather_mupexi_by_suffix.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
marker="${logdir}/submitted.${prefix}.${job_tag}.jobid"
mkdir -p "$logdir" "$repdir"

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
      echo "[skip-running] ${prefix}.${job_tag}: active job ${prev_jobid} (${st})"
      exit 0
    fi
  fi
fi

active_jobid=""
if command -v qselect >/dev/null 2>&1; then
  active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "${prefix}.${job_tag}" 2>/dev/null | head -n1 || true)"
fi
if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
  active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="${prefix}.${job_tag}" '$4==n {print $1; exit}')"
fi
if [ -n "$active_jobid" ]; then
  echo "[skip] ${prefix}.${job_tag}: scheduler already has active job ${active_jobid}"
  exit 0
fi

runscript="${logdir}/run.${job_tag}.${prefix}.sh"
patient_snippet=""
if [ "${#patient_arg[@]}" -gt 0 ]; then
  patient_snippet="$(printf ' \\\n  --patient %s' "$(printf '%q' "$patient_tag")")"
fi

cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${repo_dir}/research/gather_mupexi_by_suffix.py" \\
  --input-dir $(printf '%q' "$input_dir") \\
  --suffix $(printf '%q' "$suffix") \\
  --outfile $(printf '%q' "$outfile")${patient_snippet}

if [ ! -s $(printf '%q' "$outfile") ]; then
  echo "ERROR: gather_mupexi_by_suffix output missing/empty: ${outfile}" >&2
  exit 2
fi
SCRIPT
chmod +x "$runscript"

qsub_opts=()
[ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
[ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
qsub_opts+=(-N "${prefix}.${job_tag}")
qsub_opts+=(-o "${repdir}/${prefix}.${job_tag}.o\$PBS_JOBID")
qsub_opts+=(-e "${repdir}/${prefix}.${job_tag}.e\$PBS_JOBID")

jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
printf '%s\n' "$jobid" > "$marker"
echo "[submit] ${prefix}.${job_tag}: suffix=${suffix} jobid=${jobid}"
echo "[info] output: ${outfile}"
echo ".. logs and reports saved in ${logroot}"
