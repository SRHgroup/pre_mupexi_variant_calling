#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash splicing/run_liftover_normal_reference.sh [-c CONFIG] --input REF37.tsv.gz --out REF38.tsv.gz --chain hg19ToHg38.over.chain.gz [--engine python|ucsc] [--liftover-bin liftOver] [--from-build GRCh37] [--to-build GRCh38] [-f] [--dry-run]

Behavior:
- Submits one PBS/qsub liftover job
- Lifts the two 1-bp splice-boundary positions in a compact normal reference.
- The default Python engine parses the UCSC chain file directly and does not need a liftOver binary.
- Writes a GRCh38-compatible compact normal reference plus unmapped/summary audit files
USAGE
}

config=""
input_ref=""
out_ref=""
chain_file=""
liftover_engine=""
liftover_bin=""
from_build="GRCh37"
to_build="GRCh38"
force=0
dry_run=0
keep_cross_chrom=0
allow_inverted=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="${2:-}"; shift 2 ;;
    --input) input_ref="${2:-}"; shift 2 ;;
    --out) out_ref="${2:-}"; shift 2 ;;
    --chain) chain_file="${2:-}"; shift 2 ;;
    --engine) liftover_engine="${2:-}"; shift 2 ;;
    --liftover-bin) liftover_bin="${2:-}"; shift 2 ;;
    --from-build) from_build="${2:-}"; shift 2 ;;
    --to-build) to_build="${2:-}"; shift 2 ;;
    --keep-cross-chrom) keep_cross_chrom=1; shift ;;
    --allow-inverted) allow_inverted=1; shift ;;
    -f|--force) force=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
script_path="${repo_root}/splicing/liftover_normal_reference.py"
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

[ -n "$input_ref" ] || input_ref="${splicing_normal_junction_ref_grch37:-}"
[ -n "$out_ref" ] || out_ref="${splicing_normal_junction_ref_grch38:-}"
[ -n "$chain_file" ] || chain_file="${splicing_liftover_chain:-}"
[ -n "$liftover_engine" ] || liftover_engine="${splicing_liftover_engine:-python}"
[ -n "$liftover_bin" ] || liftover_bin="${splicing_liftover_bin:-liftOver}"

[ -n "$input_ref" ] || { echo "ERROR: --input is required unless CONFIG defines splicing_normal_junction_ref_grch37" >&2; exit 1; }
[ -n "$out_ref" ] || { echo "ERROR: --out is required unless CONFIG defines splicing_normal_junction_ref_grch38" >&2; exit 1; }
[ -n "$chain_file" ] || { echo "ERROR: --chain is required unless CONFIG defines splicing_liftover_chain" >&2; exit 1; }
[ -f "$input_ref" ] || { echo "ERROR: input reference not found: $input_ref" >&2; exit 1; }
[ -f "$chain_file" ] || { echo "ERROR: chain file not found: $chain_file" >&2; exit 1; }

if [ -f "$out_ref" ] && [ "$force" != "1" ]; then
  echo "[skip] output exists: $out_ref"
  echo "       pass -f to rebuild"
  exit 0
fi

out_dir="$(cd "$(dirname "$out_ref")" && pwd 2>/dev/null || true)"
if [ -z "$out_dir" ]; then
  out_dir="$(dirname "$out_ref")"
fi
mkdir -p "$out_dir"

logroot="${splicing_liftover_logroot:-${out_dir}/normal_reference_liftover.logs_and_reports}"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

job_name="${splicing_liftover_job_name:-splicing_liftover_ref}"
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
for module_name in \${splicing_liftover_modules:-\${splicing_reference_modules:-tools \${splicing_python_modules:-anaconda3/2025.06-1}}}; do
  module load "\$module_name"
done
splicing_python="\${splicing_python:-python3}"
printf '[liftover-normal-ref] Python: %s\\n' "\$(command -v "\$splicing_python" || printf '%s' "\$splicing_python")"
printf '[liftover-normal-ref] engine: %s\\n' $(printf '%q' "$liftover_engine")
cmd=("\$splicing_python" $(printf '%q' "$script_path") --input $(printf '%q' "$input_ref") --out $(printf '%q' "$out_ref") --chain $(printf '%q' "$chain_file") --engine $(printf '%q' "$liftover_engine") --liftover-bin $(printf '%q' "$liftover_bin") --from-build $(printf '%q' "$from_build") --to-build $(printf '%q' "$to_build"))
if [ "$keep_cross_chrom" -eq 1 ]; then
  cmd+=(--keep-cross-chrom)
fi
if [ "$allow_inverted" -eq 1 ]; then
  cmd+=(--allow-inverted)
fi
"\${cmd[@]}"
SCRIPT
chmod +x "$runscript"

if [ "$dry_run" = "1" ]; then
  echo "[dry-run] would submit ${job_name}: $runscript"
  echo "[dry-run] qsub -N ${job_name} -o ${repdir}/${job_name}.o\\\$PBS_JOBID -e ${repdir}/${job_name}.e\\\$PBS_JOBID $runscript"
  exit 0
fi

qsub_resources="nodes=${splicing_liftover_qsub_nodes:-1}:ppn=${splicing_liftover_qsub_ppn:-2},mem=${splicing_liftover_qsub_mem:-32gb},walltime=${splicing_liftover_qsub_walltime:-24:00:00}"
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
