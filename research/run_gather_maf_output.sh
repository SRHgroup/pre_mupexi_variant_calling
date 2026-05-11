#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_gather_maf_output.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [-f] [--skip-running]
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
: "${vcfdir:?CONFIG must define vcfdir}"

if [ -z "$outdir" ]; then
  if [ -n "${mupexi_outdir:-}" ]; then
    outdir="${mupexi_outdir%/}/gathered"
  else
    outdir="${vcfdir%/}/maf_output"
  fi
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
rna_tumor_sample_name="${mupexi_tumor_sample:-${rna7_signal_sample_label:-TUMOR}}"
dna_tumor_sample_name="${mupexi_dna_only_tumor_sample:-${rna_tumor_sample_name}}"
phased_ext="${rna7_phased_vcf_extension:-${phased_vcf_extension:-}}"
dna_only_phased_ext="${dna_only_phased_vcf_extension:-}"
vep_dir="${gather_maf_vep_dir:-${vep_dedup_vep_dir:-${variant_table_vep_dir:-${rna_edit_vep_dir:-${mupexi_outdir:-${datadir%/}/mupexi2}}}}}"

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

patient_has_rna_sample() {
  local patient="$1"
  local sid base label
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
    case "${sid,,}" in \#sample|sample|sample_id|patient|patient_id) continue ;; esac
    base="$(sample_base_name "$sid")"
    [ "$base" = "$patient" ] || continue
    label="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $4}')"
    if [ "$label" = "${rna_tumor_label:-RNA_TUMOR}" ] || [ "$label" = "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}" ]; then
      return 0
    fi
    if [[ "$sid" == *"_${rna_tumor_label:-RNA_TUMOR}" ]] || [[ "$sid" == *"_${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}" ]] || [[ "$sid" == *"_RNA_TUMOR" ]] || [[ "$sid" == *"_RNA_TUMOUR" ]]; then
      return 0
    fi
  done < "$samples"
  return 1
}

vep_inputs=()
vcf_inputs=()
tumor_samples=()
seen=""
patient_tag="cohort"
while IFS= read -r line; do
  [ -n "$line" ] || continue
  [[ "$line" =~ ^[[:space:]]*# ]] && continue
  sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
  case "${sid,,}" in \#sample|sample|sample_id|patient|patient_id) continue ;; esac
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

  vep="${vep_dir}/${patient}_vep.vep"
  if [ ! -f "$vep" ] && [ -f "${vep}.gz" ]; then
    vep="${vep}.gz"
  fi
  if [ ! -f "$vep" ]; then
    echo "[skip] ${patient}: missing VEP file: ${vep_dir}/${patient}_vep.vep(.gz)"
    continue
  fi

  vcf=""
  tumor_sample="$dna_tumor_sample_name"
  if patient_has_rna_sample "$patient"; then
    vcf="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}/${patient}_${phased_ext}"
    tumor_sample="$rna_tumor_sample_name"
  else
    vcf="${vcfdir}/${patient}_${out_dna_label}_vs_${patient}_${out_normal_label}/${patient}_${dna_only_phased_ext}"
  fi
  if [ ! -f "$vcf" ]; then
    echo "[skip] ${patient}: missing phased VCF: $vcf"
    continue
  fi

  vep_inputs+=(--vep-input "${patient}=${vep}")
  vcf_inputs+=(--vcf-input "${patient}=${vcf}")
  tumor_samples+=(--tumor-sample "${patient}=${tumor_sample}")
done < "$samples"

if [ "${#vep_inputs[@]}" -eq 0 ]; then
  echo "ERROR: no VEP/phased VCF input pairs found" >&2
  exit 1
fi

outfile="${outdir}/${patient_tag}.maf_like.tsv"
if [ "$force" != "1" ] && [ -s "$outfile" ]; then
  echo "[skip] output already exists: $outfile (use -f to overwrite)"
  exit 0
fi

logroot="${outdir}/gather_maf_output.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_gather_maf_output"
marker="${logdir}/submitted.${prefix}.${patient_tag}.jobid"

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
      echo "[skip-running] ${prefix}.${patient_tag}: active job ${prev_jobid} (${st})"
      exit 0
    fi
  fi
fi

active_jobid=""
if command -v qselect >/dev/null 2>&1; then
  active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "${prefix}.${patient_tag}" 2>/dev/null | head -n1 || true)"
fi
if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
  active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="${prefix}.${patient_tag}" '$4==n {print $1; exit}')"
fi
if [ -n "$active_jobid" ]; then
  echo "[skip] ${prefix}.${patient_tag}: scheduler already has active job ${active_jobid}"
  exit 0
fi

runscript="${logdir}/run.${patient_tag}.${prefix}.sh"
apply_vep_inputs=()
for item in "${vep_inputs[@]}"; do
  apply_vep_inputs+=("$(printf '%q' "$item")")
done
apply_vcf_inputs=()
for item in "${vcf_inputs[@]}"; do
  apply_vcf_inputs+=("$(printf '%q' "$item")")
done
apply_tumor_samples=()
for item in "${tumor_samples[@]}"; do
  apply_tumor_samples+=("$(printf '%q' "$item")")
done

cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${repo_dir}/research/gather_maf_output.py" \\
  ${apply_vep_inputs[*]} \\
  ${apply_vcf_inputs[*]} \\
  ${apply_tumor_samples[*]} \\
  --tumor-label "$(printf '%q' "${rna7_signal_sample_label:-TUMOR}")" \\
  --tumor-label "$(printf '%q' "$out_dna_label")" \\
  --normal-label "$(printf '%q' "$out_normal_label")" \\
  --outfile "$(printf '%q' "$outfile")"

if [ ! -s "$(printf '%q' "$outfile")" ]; then
  echo "ERROR: gather_maf_output output missing/empty: $(printf '%q' "$outfile")" >&2
  exit 2
fi
SCRIPT
chmod +x "$runscript"

qsub_opts=()
[ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
[ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
qsub_opts+=(-N "${prefix}.${patient_tag}")
qsub_opts+=(-o "${repdir}/${prefix}.${patient_tag}.o\$PBS_JOBID")
qsub_opts+=(-e "${repdir}/${prefix}.${patient_tag}.e\$PBS_JOBID")

jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
printf '%s\n' "$jobid" > "$marker"
echo "[submit] ${prefix}.${patient_tag}: jobid=${jobid}"
echo "[info] gather_maf_output file: ${outfile}"
echo ".. logs and reports saved in ${logroot}"
