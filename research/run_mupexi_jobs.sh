#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_mupexi_jobs.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [--run-fusions] [-f] [--skip-running]
USAGE
}

config=""
sample=""
outdir=""
run_fusions=0
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --run-fusions) run_fusions=1; shift ;;
    -f|--force) force=1; shift ;;
    --skip-running) skip_running=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[ -n "$config" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${mupexi2_repo:?CONFIG must define mupexi2_repo}"
: "${mupexi_netmhc_config:?CONFIG must define mupexi_netmhc_config}"

if [ -z "$outdir" ]; then
  outdir="${mupexi_outdir:-${vcfdir}/mupexi2}"
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
tumor_sample_name="${mupexi_tumor_sample:-TUMOR}"
normal_sample_name="${mupexi_normal_sample:-DNA_NORMAL}"
peptide_lengths="${mupexi_peptide_lengths:-9-11}"
parallel_k="${mupexi_parallel_k:-true}"
enable_germlines="${mupexi_enable_germlines:-true}"
enable_superpeptides="${mupexi_enable_superpeptides:-true}"
enable_rna_edit="${mupexi_enable_rna_edit:-true}"
phased_ext="${rna7_phased_vcf_extension:-${phased_vcf_extension:-}}"
[ -n "$phased_ext" ] || { echo "ERROR: missing phased VCF extension in CONFIG (rna7_phased_vcf_extension/phased_vcf_extension)" >&2; exit 1; }

resolve_patient_placeholder() {
  local template="$1"
  local patient="$2"
  printf '%s\n' "${template//\{patient\}/$patient}"
}

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

lookup_map_value() {
  local map_file="$1"
  local patient="$2"
  [ -f "$map_file" ] || return 1
  awk -F'\t' -v p="$patient" '$1==p {print $2; exit}' "$map_file"
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

logroot="${outdir}/mupexi.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

declare -A seen=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac
  sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
  patient="$(sample_base_name "$sid")"
  [ -n "$patient" ] || continue
  [ -n "${seen[$patient]:-}" ] && continue
  seen["$patient"]=1
  if [ -n "$sample" ] && [ "$sample" != "$sid" ] && [ "$sample" != "$patient" ]; then
    continue
  fi

  vcf="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}/${patient}_${phased_ext}"
  if [ ! -f "$vcf" ]; then
    echo "[skip] ${patient}: missing phased VCF: $vcf"
    continue
  fi

  expr=""
  if [ -n "${mupexi_expression_map_tsv:-}" ]; then
    expr="$(lookup_map_value "$mupexi_expression_map_tsv" "$patient" || true)"
  elif [ -n "${mupexi_expression_tsv_template:-}" ]; then
    expr="$(resolve_patient_placeholder "$mupexi_expression_tsv_template" "$patient")"
  fi
  if [ -z "$expr" ] || [ ! -f "$expr" ]; then
    echo "[skip] ${patient}: missing expression file (set mupexi_expression_tsv_template or mupexi_expression_map_tsv)"
    continue
  fi

  hla=""
  if [ -n "${mupexi_hla_map_tsv:-}" ]; then
    hla="$(lookup_map_value "$mupexi_hla_map_tsv" "$patient" || true)"
  fi
  if [ -z "$hla" ]; then
    hla="${mupexi_hla:-}"
  fi
  if [ -z "$hla" ]; then
    echo "[skip] ${patient}: missing HLA string (set mupexi_hla or mupexi_hla_map_tsv)"
    continue
  fi

  fusion_arg=""
  if [ "$run_fusions" = "1" ]; then
    fusion_path=""
    if [ -n "${mupexi_fusion_arriba_map_tsv:-}" ]; then
      fusion_path="$(lookup_map_value "$mupexi_fusion_arriba_map_tsv" "$patient" || true)"
    elif [ -n "${mupexi_fusion_arriba_template:-}" ]; then
      fusion_path="$(resolve_patient_placeholder "$mupexi_fusion_arriba_template" "$patient")"
    fi
    if [ -n "$fusion_path" ] && [ -f "$fusion_path" ]; then
      fusion_arg="--fusion ${fusion_path}"
    else
      echo "[warn] ${patient}: --run-fusions requested but fusion_arriba file missing; running without -z"
    fi
  fi

  # Heuristic skip if patient-prefixed outputs already exist.
  if [ "$force" != "1" ] && find "$outdir" -maxdepth 1 \( -type f -o -type d \) -name "${patient}*" | grep -q .; then
    echo "[skip] ${patient}: mupexi output(s) already exist in ${outdir} (use -f to overwrite)"
    continue
  fi

  prefix="mupexi2"
  marker="${logdir}/submitted.${prefix}.${patient}.jobid"
  if [ "$skip_running" = "1" ] && [ -f "$marker" ]; then
    prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ]; then
      st="$(pbs_state_for_jobid "$prev_jobid")"
      if [ "$st" = "RUNNING" ] || [ "$st" = "QUEUED" ]; then
        echo "[skip-running] ${prefix}.${patient}: active job ${prev_jobid} (${st})"
        continue
      fi
    fi
  fi

  runscript="${logdir}/run.${patient}.${prefix}.sh"
  cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${mupexi_modules:-tools ngs anaconda3/2025.06-1 netmhcpan/4.0a perl/5.36.1 ensembl-tools/90}

export PYTHONPATH="${mupexi2_repo}/src:\${PYTHONPATH:-}"

python3 -m mupexi2.cli \\
  -v "${vcf}" \\
  --vcf-type merged \\
  --germlines "${enable_germlines}" \\
  --superpeptides "${enable_superpeptides}" \\
  --rna-edit "${enable_rna_edit}" \\
  --tumor-sample "${tumor_sample_name}" \\
  --normal-sample "${normal_sample_name}" \\
  --parallel-k "${parallel_k}" \\
  -l "${peptide_lengths}" \\
  -a "${hla}" \\
  -t -f -n \\
  -c "${mupexi_netmhc_config}" \\
  -p "${patient}" \\
  -d "${outdir}" \\
  -e "${expr}" \\
  ${fusion_arg}
SCRIPT
  chmod +x "$runscript"

  qsub_opts=()
  [ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
  [ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
  qsub_opts+=(-N "${prefix}.${patient}")
  qsub_opts+=(-o "${repdir}/${prefix}.${patient}.o\$PBS_JOBID")
  qsub_opts+=(-e "${repdir}/${prefix}.${patient}.e\$PBS_JOBID")

  jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
  printf '%s\n' "$jobid" > "$marker"
  echo "[submit] ${prefix}.${patient}: jobid=${jobid}"
done < "$samples"

echo ".. logs and reports saved in ${logroot}"
