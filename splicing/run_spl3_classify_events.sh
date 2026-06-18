#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash splicing/run_spl3_classify_events.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [--root SPLICING_DIR] [--gtf GTF] [--outdir DIR] [-f] [--dry-run] [--include-noncanonical]

Behavior:
- Submits a PBS/qsub job per patient/sample
- The qsub job reads spl2 *.spl2.novel_junctions.tsv files from the splicing root
- The qsub job classifies each junction against GTF transcript exon/intron structure
- Writes *.spl3.event_annotated.tsv under --outdir, CONFIG:splicing_outdir,
  or ${datadir}/splicing when datadir is defined
USAGE
}

config=""
sample=""
root_override=""
gtf_override=""
outdir_override=""
force=0
dry_run=0
include_noncanonical=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="${2:-}"; shift 2 ;;
    -s|--sample) sample="${2:-}"; shift 2 ;;
    --root) root_override="${2:-}"; shift 2 ;;
    --gtf) gtf_override="${2:-}"; shift 2 ;;
    --outdir) outdir_override="${2:-}"; shift 2 ;;
    -f|--force) force=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    --include-noncanonical) include_noncanonical=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
script_path="${repo_root}/splicing/classify_junction_events.py"
pipeline_defaults="${PIPELINE_DEFAULTS:-${repo_root}/pipeline_defaults/toolchain.defaults.sh}"

[ -n "$config" ] || { usage >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

if [ -n "${pipeline_defaults:-}" ] && [ -f "$pipeline_defaults" ]; then
  # shellcheck disable=SC1090
  source "$pipeline_defaults"
fi

# shellcheck disable=SC1090
source "$config"

if [ -n "$root_override" ]; then
  spl3_root="$root_override"
elif [ -n "${splicing_outdir:-}" ]; then
  spl3_root="$splicing_outdir"
elif [ -n "${datadir:-}" ]; then
  spl3_root="${datadir%/}/splicing"
else
  echo "ERROR: CONFIG must define splicing_outdir or datadir; or pass --root" >&2
  exit 1
fi

if [ -n "$gtf_override" ]; then
  gtf_path="$gtf_override"
elif [ -n "${GTF:-}" ]; then
  gtf_path="$GTF"
else
  echo "ERROR: CONFIG must define GTF, or pass --gtf" >&2
  exit 1
fi

if [ -n "$outdir_override" ]; then
  spl3_outdir="$outdir_override"
elif [ -n "${splicing_outdir:-}" ]; then
  spl3_outdir="$splicing_outdir"
elif [ -n "${datadir:-}" ]; then
  spl3_outdir="${datadir%/}/splicing"
else
  spl3_outdir="$spl3_root"
fi

[ -d "$spl3_root" ] || { echo "ERROR: spl2/splicing root not found: $spl3_root" >&2; exit 1; }
[ -f "$gtf_path" ] || { echo "ERROR: missing GTF: $gtf_path" >&2; exit 1; }

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

is_rna_sample_id() {
  local sid="$1"
  case "$sid" in
    *"_${rna_tumor_label:-RNA_TUMOR}"|*"_${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"|*_RNA_TUMOR|*_RNA_TUMOUR) return 0 ;;
    *) return 1 ;;
  esac
}

discover_targets() {
  if [ -n "$sample" ]; then
    printf '%s\n' "$sample"
    return 0
  fi

  : "${samples:?CONFIG must define samples when no sample is requested}"
  [ -f "$samples" ] || { echo "ERROR: samples file not found: $samples" >&2; exit 1; }

  local seen="" line sid patient
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    case "$line" in [[:space:]]*'#'*) continue ;; esac
    sid="$(printf '%s\n' "$line" | awk -F'[,	 ]+' '{print $1}')"
    [ -n "$sid" ] || continue
    is_rna_sample_id "$sid" || continue
    patient="$(sample_base_name "$sid")"
    [ -n "$patient" ] || continue
    if printf '%s\n' "$seen" | grep -Fxq "$patient"; then
      continue
    fi
    seen="${seen}
${patient}"
    printf '%s\n' "$patient"
  done < "$samples"
}

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
  local job_name="$1"
  local active_jobid=""
  if command -v qselect >/dev/null 2>&1; then
    active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "$job_name" 2>/dev/null | head -n1 || true)"
  fi
  if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
    active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="$job_name" '$4==n {print $1; exit}')"
  fi
  printf '%s\n' "$active_jobid"
}

submit_target() {
  local target="$1"
  local tag="${target//[^A-Za-z0-9_.-]/_}"
  local prefix="splicing_spl3"
  local job_name="${prefix}.${tag}"
  local logroot="${splicing_logroot:-${spl3_outdir:-${spl3_root%/}}/${prefix}.logs_and_reports}"
  local logdir="${logroot}/logs"
  local repdir="${logroot}/reports"
  local marker="${logdir}/submitted.${job_name}.jobid"
  local runscript="${logdir}/run.${tag}.${prefix}.sh"
  mkdir -p "$logdir" "$repdir"

  if [ "$dry_run" != "1" ] && [ -f "$marker" ]; then
    local prev_jobid st
    prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ]; then
      st="$(pbs_state_for_jobid "$prev_jobid")"
      if [ "$st" = "RUNNING" ] || [ "$st" = "QUEUED" ]; then
        echo "[skip] ${job_name}: active job ${prev_jobid} (${st})"
        return 0
      fi
    fi
  fi

  if [ "$dry_run" != "1" ]; then
    local active_jobid
    active_jobid="$(active_job_for_name "$job_name")"
    if [ -n "$active_jobid" ]; then
      echo "[skip] ${job_name}: scheduler already has active job ${active_jobid}"
      return 0
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
# shellcheck disable=SC1090
source $(printf '%q' "$config")
if [ -n "\${splicing_module_prereq:-tools}" ]; then
  module load \${splicing_module_prereq:-tools}
fi
module load \${splicing_python_modules:-anaconda3/2025.06-1}
splicing_python="\${splicing_python:-python3}"
printf '[spl3] spl2 root: %s\\n' $(printf '%q' "$spl3_root")
printf '[spl3] GTF: %s\\n' $(printf '%q' "$gtf_path")
spl3_outdir=$(printf '%q' "$spl3_outdir")
if [ -n "\$spl3_outdir" ]; then
  mkdir -p "\$spl3_outdir"
  printf '[spl3] output root: %s\\n' "\$spl3_outdir"
fi
printf '[spl3] Python: %s\\n' "\$(command -v "\$splicing_python" || printf '%s' "\$splicing_python")"
cmd=("\$splicing_python" $(printf '%q' "$script_path") --root $(printf '%q' "$spl3_root") --gtf $(printf '%q' "$gtf_path") --sample-filter $(printf '%q' "$target"))
if [ -n "\$spl3_outdir" ]; then
  cmd+=(--out-dir "\$spl3_outdir")
fi
if [ "$force" -eq 1 ]; then
  cmd+=(--force)
fi
if [ "$include_noncanonical" -eq 1 ]; then
  cmd+=(--include-noncanonical)
fi
"\${cmd[@]}"
SCRIPT
  chmod +x "$runscript"

  if [ "$dry_run" = "1" ]; then
    echo "[dry-run] would submit ${job_name}: $runscript"
    echo "[dry-run] qsub -N ${job_name} -o ${repdir}/${job_name}.o\\\$PBS_JOBID -e ${repdir}/${job_name}.e\\\$PBS_JOBID $runscript"
    return 0
  fi

  local qsub_resources="nodes=${splicing_spl3_qsub_nodes:-${splicing_qsub_nodes:-1}}:ppn=${splicing_spl3_qsub_ppn:-${splicing_qsub_ppn:-2}},mem=${splicing_spl3_qsub_mem:-${splicing_qsub_mem:-16gb}},walltime=${splicing_spl3_qsub_walltime:-${splicing_qsub_walltime:-06:00:00}}"
  local qsub_opts=()
  [ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
  [ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
  qsub_opts+=(-l "$qsub_resources")
  qsub_opts+=(-N "$job_name")
  qsub_opts+=(-o "${repdir}/${job_name}.o\$PBS_JOBID")
  qsub_opts+=(-e "${repdir}/${job_name}.e\$PBS_JOBID")

  local jobid
  jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
  printf '%s\n' "$jobid" > "$marker"
  echo "[submit] ${job_name}: jobid=${jobid}"
  echo ".. logs and reports saved in ${logroot}"
}

submitted=0
while IFS= read -r target; do
  [ -n "$target" ] || continue
  submit_target "$target"
  submitted=$((submitted + 1))
done < <(discover_targets)

if [ "$submitted" -eq 0 ]; then
  echo "ERROR: no RNA patient/sample targets found for splicing spl3" >&2
  exit 1
fi
