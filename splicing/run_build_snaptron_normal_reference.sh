#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash splicing/run_build_snaptron_normal_reference.sh [-c CONFIG] --snaptron-dir DIR [--out TSV.GZ] [--coordinate-mode star-intron|boundary] [--total-samples N] [--min-sample-count N] [--min-total-reads N] [--min-prevalence X] [--canonical-only] [--drop-unknown-strand] [--include-tissue-summary] [--sample-filter-column COLUMN=VALUE] [-f] [--dry-run]

Behavior:
- Submits one PBS/qsub reference-build job
- Input DIR must contain junctions.bgz, samples.tsv, and samples.fields.tsv
- Output is a compact normal junction TSV usable by splicing spl3.5
USAGE
}

config=""
snaptron_dir=""
out_path=""
coordinate_mode="star-intron"
total_samples=""
min_sample_count="1"
min_total_reads="1"
min_prevalence="0"
canonical_only=0
drop_unknown_strand=0
include_tissue_summary=0
force=0
dry_run=0
sample_filter_columns=()

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="${2:-}"; shift 2 ;;
    --snaptron-dir) snaptron_dir="${2:-}"; shift 2 ;;
    --out) out_path="${2:-}"; shift 2 ;;
    --coordinate-mode) coordinate_mode="${2:-}"; shift 2 ;;
    --total-samples) total_samples="${2:-}"; shift 2 ;;
    --min-sample-count) min_sample_count="${2:-}"; shift 2 ;;
    --min-total-reads) min_total_reads="${2:-}"; shift 2 ;;
    --min-prevalence) min_prevalence="${2:-}"; shift 2 ;;
    --canonical-only) canonical_only=1; shift ;;
    --drop-unknown-strand) drop_unknown_strand=1; shift ;;
    --include-tissue-summary) include_tissue_summary=1; shift ;;
    --sample-filter-column) sample_filter_columns+=("${2:-}"); shift 2 ;;
    -f|--force) force=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
script_path="${repo_root}/splicing/build_snaptron_normal_reference.py"
pipeline_defaults="${PIPELINE_DEFAULTS:-${repo_root}/pipeline_defaults/toolchain.defaults.sh}"

if [ -n "${pipeline_defaults:-}" ] && [ -f "$pipeline_defaults" ]; then
  # shellcheck disable=SC1090
  source "$pipeline_defaults"
fi

if [ -n "$config" ]; then
  [ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$config"
fi

[ -n "$snaptron_dir" ] || { echo "ERROR: --snaptron-dir is required" >&2; usage >&2; exit 1; }
[ -d "$snaptron_dir" ] || { echo "ERROR: Snaptron folder not found: $snaptron_dir" >&2; exit 1; }
[ -f "${snaptron_dir%/}/junctions.bgz" ] || { echo "ERROR: missing ${snaptron_dir%/}/junctions.bgz" >&2; exit 1; }
[ -f "${snaptron_dir%/}/samples.tsv" ] || { echo "ERROR: missing ${snaptron_dir%/}/samples.tsv" >&2; exit 1; }
[ -f "${snaptron_dir%/}/samples.fields.tsv" ] || { echo "ERROR: missing ${snaptron_dir%/}/samples.fields.tsv" >&2; exit 1; }

if [ -z "$out_path" ]; then
  out_path="${splicing_normal_junction_ref:-}"
fi
[ -n "$out_path" ] || { echo "ERROR: --out is required unless CONFIG defines splicing_normal_junction_ref" >&2; exit 1; }

if [ -z "$total_samples" ]; then
  total_samples="${splicing_normal_total_samples:-}"
fi

if [ -f "$out_path" ] && [ "$force" != "1" ]; then
  echo "[skip] output exists: $out_path"
  echo "       pass -f to rebuild"
  exit 0
fi

out_dir="$(cd "$(dirname "$out_path")" && pwd 2>/dev/null || true)"
if [ -z "$out_dir" ]; then
  out_dir="$(dirname "$out_path")"
fi
mkdir -p "$out_dir"

logroot="${splicing_reference_logroot:-${out_dir}/snaptron_normal_reference.logs_and_reports}"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

job_name="${splicing_reference_job_name:-splicing_snaptron_ref}"
runscript="${logdir}/run.${job_name}.sh"
marker="${logdir}/submitted.${job_name}.jobid"

pbs_state_for_jobid() {
  local jid="$1"
  local state
  state="$(qstat -f "$jid" 2>/dev/null | awk -F' = ' '/job_state =/{print $2; exit}' || true)"
  case "$state" in
    R|E) printf '%s\n' "RUNNING" ;;
    Q|H|W|T|S) printf '%s\n' "QUEUED" ;;
    *) printf '%s\n' "" ;;
  esac
}

active_job_for_name() {
  local active_jobid=""
  if command -v qselect >/dev/null 2>&1; then
    active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "$job_name" 2>/dev/null | head -n1 || true)"
  fi
  if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
    active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="$job_name" '$4==n {print $1; exit}')"
  fi
  printf '%s\n' "$active_jobid"
}

if [ "$dry_run" != "1" ] && [ -f "$marker" ]; then
  prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
  if [ -n "$prev_jobid" ]; then
    st="$(pbs_state_for_jobid "$prev_jobid")"
    if [ "$st" = "RUNNING" ] || [ "$st" = "QUEUED" ]; then
      echo "[skip] ${job_name}: active job ${prev_jobid} (${st})"
      exit 0
    fi
  fi
fi

if [ "$dry_run" != "1" ]; then
  active_jobid="$(active_job_for_name)"
  if [ -n "$active_jobid" ]; then
    echo "[skip] ${job_name}: scheduler already has active job ${active_jobid}"
    exit 0
  fi
fi

sample_filter_snippet=""
for sample_filter in "${sample_filter_columns[@]}"; do
  sample_filter_snippet+="cmd+=(--sample-filter-column $(printf '%q' "$sample_filter"))"$'\n'
done

cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
export PIPELINE_DEFAULTS=$(printf '%q' "$pipeline_defaults")
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
if [ -n $(printf '%q' "$config") ]; then
  # shellcheck disable=SC1090
  source $(printf '%q' "$config")
fi
for module_name in \${splicing_reference_modules:-tools \${splicing_python_modules:-anaconda3/2025.06-1}}; do
  module load "\$module_name"
done
splicing_python="\${splicing_python:-python3}"
printf '[snaptron-ref] Python: %s\\n' "\$(command -v "\$splicing_python" || printf '%s' "\$splicing_python")"
cmd=("\$splicing_python" $(printf '%q' "$script_path") --snaptron-dir $(printf '%q' "$snaptron_dir") --out $(printf '%q' "$out_path") --coordinate-mode $(printf '%q' "$coordinate_mode") --min-sample-count $(printf '%q' "$min_sample_count") --min-total-reads $(printf '%q' "$min_total_reads") --min-prevalence $(printf '%q' "$min_prevalence"))
if [ -n $(printf '%q' "$total_samples") ]; then
  cmd+=(--total-samples $(printf '%q' "$total_samples"))
fi
if [ "$canonical_only" -eq 1 ]; then
  cmd+=(--canonical-only)
fi
if [ "$drop_unknown_strand" -eq 1 ]; then
  cmd+=(--drop-unknown-strand)
fi
if [ "$include_tissue_summary" -eq 1 ]; then
  cmd+=(--include-tissue-summary)
fi
${sample_filter_snippet}"\${cmd[@]}"
SCRIPT
chmod +x "$runscript"

if [ "$dry_run" = "1" ]; then
  echo "[dry-run] would submit ${job_name}: $runscript"
  echo "[dry-run] qsub -N ${job_name} -o ${repdir}/${job_name}.o\\\$PBS_JOBID -e ${repdir}/${job_name}.e\\\$PBS_JOBID $runscript"
  exit 0
fi

qsub_resources="nodes=${splicing_reference_qsub_nodes:-1}:ppn=${splicing_reference_qsub_ppn:-2},mem=${splicing_reference_qsub_mem:-32gb},walltime=${splicing_reference_qsub_walltime:-24:00:00}"
qsub_opts=()
[ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
[ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
qsub_opts+=(-l "$qsub_resources")
qsub_opts+=(-N "$job_name")
qsub_opts+=(-o "${repdir}/${job_name}.o\$PBS_JOBID")
qsub_opts+=(-e "${repdir}/${job_name}.e\$PBS_JOBID")

jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
printf '%s\n' "$jobid" > "$marker"
echo "[submit] ${job_name}: jobid=${jobid}"
echo ".. logs and reports saved in ${logroot}"
