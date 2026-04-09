#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash bin/run_cleanup_pre_mupexi_jobs.sh -c CONFIG [-s PATIENT] [--execute] [--threads N] [--nodes N] [--ppn N] [--mem SIZE] [--walltime HH:MM:SS] [-f] [--skip-running]

Behavior:
  - without --execute: run local dry-run planning via cleanup_pre_mupexi_outputs.sh
  - with --execute: submit one PBS job per ready patient
USAGE
}

config=""
sample=""
execute=0
force=0
skip_running=0
threads=""
cli_nodes=""
cli_ppn=""
cli_mem=""
cli_walltime=""

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    --execute) execute=1; shift ;;
    --threads) threads="$2"; shift 2 ;;
    --nodes) cli_nodes="$2"; shift 2 ;;
    --ppn) cli_ppn="$2"; shift 2 ;;
    --mem) cli_mem="$2"; shift 2 ;;
    --walltime) cli_walltime="$2"; shift 2 ;;
    -f|--force) force=1; shift ;;
    --skip-running) skip_running=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$config" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "$PIPELINE_DEFAULTS"
fi

# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
: "${bamdir:?CONFIG must define bamdir}"
[ -f "$samples" ] || { echo "ERROR: samples file not found: $samples" >&2; exit 1; }

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
worker="${repo_dir}/bin/cleanup_pre_mupexi_outputs.sh"
[ -f "$worker" ] || { echo "ERROR: cleanup worker not found: $worker" >&2; exit 1; }

prefix="cleanup_pre_mupexi"
logroot="${bamdir}/${prefix}.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

q_nodes="${cli_nodes:-${cleanup_qsub_nodes:-1}}"
q_ppn="${cli_ppn:-${cleanup_qsub_ppn:-8}}"
q_mem="${cli_mem:-${cleanup_qsub_mem:-32gb}}"
q_walltime="${cli_walltime:-${cleanup_qsub_walltime:-24:00:00}}"
threads="${threads:-$q_ppn}"

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

if [ "$execute" -ne 1 ]; then
  if [ -n "$sample" ]; then
    exec bash "$worker" -c "$config" -s "$sample" --threads "$threads"
  fi
  exec bash "$worker" -c "$config" --threads "$threads"
fi

submitted=0
ready=0
skipped=0
seen_patients=""

while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac
  sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
  patient="$(sample_base_name "$sid")"
  [ -n "$patient" ] || continue
  if printf '%s\n' "$seen_patients" | grep -Fxq "$patient"; then
    continue
  fi
  seen_patients="${seen_patients}
${patient}"
  if [ -n "$sample" ] && [ "$sample" != "$sid" ] && [ "$sample" != "$patient" ]; then
    continue
  fi

  done_marker="${logdir}/done.${prefix}.${patient}.ok"
  job_marker="${logdir}/submitted.${prefix}.${patient}.jobid"

  if [ "$force" -ne 1 ] && [ -f "$done_marker" ]; then
    echo "[skip] ${prefix}.${patient}: already marked done ($done_marker)"
    skipped=$((skipped + 1))
    continue
  fi

  if [ "$skip_running" = "1" ] && [ -f "$job_marker" ]; then
    prev_jobid="$(head -n1 "$job_marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ]; then
      st="$(pbs_state_for_jobid "$prev_jobid")"
      if [ "$st" = "RUNNING" ] || [ "$st" = "QUEUED" ]; then
        echo "[skip-running] ${prefix}.${patient}: active job ${prev_jobid} (${st})"
        skipped=$((skipped + 1))
        continue
      fi
    fi
  fi

  active_jobid=""
  if command -v qselect >/dev/null 2>&1; then
    active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "${prefix}.${patient}" 2>/dev/null | head -n1 || true)"
  fi
  if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
    active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="${prefix}.${patient}" '$4==n {print $1; exit}')"
  fi
  if [ -n "$active_jobid" ]; then
    echo "[skip] ${prefix}.${patient}: scheduler already has active job ${active_jobid}"
    skipped=$((skipped + 1))
    continue
  fi

  tmp="$(mktemp)"
  if bash "$worker" -c "$config" -s "$patient" --threads "$threads" >"$tmp" 2>&1; then
    ready=$((ready + 1))
  else
    echo "[skip] ${prefix}.${patient}: cleanup precheck failed"
    cat "$tmp"
    rm -f "$tmp"
    skipped=$((skipped + 1))
    continue
  fi
  delete_count="$(awk -F': ' '/^  delete_count:/{print $2; exit}' "$tmp")"
  rm -f "$tmp"

  runscript="${logdir}/run.${patient}.${prefix}.sh"
  cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${modules_rna_bamfix:-tools htslib/1.23 samtools/1.23}
bash "$(printf '%q' "$worker")" -c "$(printf '%q' "$config")" -s "$(printf '%q' "$patient")" --execute --threads "$(printf '%q' "$threads")"
touch "$(printf '%q' "$done_marker")"
SCRIPT
  chmod +x "$runscript"

  qsub_opts=()
  [ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
  [ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
  [ -n "${QSUB_DEPEND:-}" ] && qsub_opts+=(-W "depend=${QSUB_DEPEND}")
  qsub_opts+=(-d "$(pwd)")
  qsub_opts+=(-l "nodes=${q_nodes}:ppn=${q_ppn},mem=${q_mem},walltime=${q_walltime}")
  qsub_opts+=(-r y)
  qsub_opts+=(-N "${prefix}.${patient}")
  qsub_opts+=(-o "${repdir}/${prefix}.${patient}.o\$PBS_JOBID")
  qsub_opts+=(-e "${repdir}/${prefix}.${patient}.e\$PBS_JOBID")

  jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
  printf '%s\n' "$jobid" > "$job_marker"
  echo "[submit] ${prefix}.${patient}: jobid=${jobid} delete_count=${delete_count:-unknown}"
  submitted=$((submitted + 1))
done < "$samples"

echo "[summary] ${prefix}: ready=${ready} submitted=${submitted} skipped=${skipped}"
echo ".. logs and reports saved in ${logroot}"
