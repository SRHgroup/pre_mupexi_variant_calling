#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_vep_dedup_job.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [-f] [--skip-running]
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
  # shellcheck disable=SC1090
  source "$PIPELINE_DEFAULTS"
fi

# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"

if [ -z "$outdir" ]; then
  outdir="${vcfdir%/}/vep_dedup"
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
phased_ext="${vep_dedup_vcf_extension:-${rna7_phased_vcf_extension:-${phased_vcf_extension:-}}}"
vep_dir="${vep_dedup_vep_dir:-${variant_table_vep_dir:-${rna_edit_vep_dir:-${mupexi_outdir:-${datadir%/}/mupexi2}}}}"
tumor_col_label="${rna7_signal_sample_label:-TUMOR}"

[ -n "$phased_ext" ] || { echo "ERROR: phased VCF extension not found in CONFIG" >&2; exit 1; }
[ -n "$vep_dir" ] || { echo "ERROR: VEP dir not found in CONFIG" >&2; exit 1; }

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

seen_patients=""
vep_inputs=()
vcf_inputs=()
patients=()
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

  vcf="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}/${patient}_${phased_ext}"
  vep="$(printf '%s' "${vep_dir}/${patient}_vep.vep")"
  if [ ! -f "$vep" ] && [ -f "${vep}.gz" ]; then
    vep="${vep}.gz"
  fi

  if [ ! -f "$vcf" ]; then
    echo "[skip] ${patient}: missing phased VCF: $vcf"
    continue
  fi
  if [ ! -f "$vep" ]; then
    echo "[skip] ${patient}: missing VEP file: ${vep_dir}/${patient}_vep.vep(.gz)"
    continue
  fi

  vep_inputs+=("--vep-input" "${patient}=${vep}")
  vcf_inputs+=("--vcf-input" "${patient}=${vcf}")
  patients+=("$patient")
done < "$samples"

if [ "${#patients[@]}" -eq 0 ]; then
  echo "ERROR: no patients with both VEP and phased VCF inputs found" >&2
  exit 1
fi

cohort_dedup="${outdir}/cohort.vep.dedup.tsv"
tag="cohort"
expected_output="$cohort_dedup"
if [ -n "$sample" ]; then
  tag="${patients[0]}"
  expected_output="${outdir}/${tag}.vep.dedup.tsv"
fi

if [ "$force" != "1" ] && [ -s "$expected_output" ]; then
  echo "[skip] output already exists: ${expected_output} (use -f to overwrite)"
  exit 0
fi

logroot="${outdir}/vep_dedup.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_vep_dedup"
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
apply_vep_inputs=()
for item in "${vep_inputs[@]}"; do
  apply_vep_inputs+=("$(printf '%q' "$item")")
done
apply_vcf_inputs=()
for item in "${vcf_inputs[@]}"; do
  apply_vcf_inputs+=("$(printf '%q' "$item")")
done

cohort_args=""
if [ -z "$sample" ]; then
  cohort_args="--cohort-dedup-outfile $(printf '%q' "$cohort_dedup") --cohort-counts-outfile $(printf '%q' "${outdir}/cohort.vep.dedup.counts.tsv") --cohort-stats-outfile $(printf '%q' "${outdir}/cohort.vep.dedup.stats.tsv")"
fi

cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${repo_dir}/research/annotate_dedup_vep.py" \\
  ${apply_vep_inputs[*]} \\
  ${apply_vcf_inputs[*]} \\
  --outdir "$(printf '%q' "$outdir")" \\
  --tumor-sample "$(printf '%q' "$tumor_col_label")" \\
  --tumor-label "$(printf '%q' "$tumor_col_label")" \\
  --tumor-label "$(printf '%q' "$out_dna_label")" \\
  --normal-label "$(printf '%q' "$out_normal_label")" \\
  ${cohort_args}

if [ ! -s "$(printf '%q' "$expected_output")" ]; then
  echo "ERROR: expected VEP dedup output missing/empty: $(printf '%q' "$expected_output")" >&2
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
echo "[info] VEP dedup output dir: ${outdir}"
echo ".. logs and reports saved in ${logroot}"
