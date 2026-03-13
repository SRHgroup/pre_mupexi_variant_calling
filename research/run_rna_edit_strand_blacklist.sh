#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_rna_edit_strand_blacklist.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [--protocol NAME] [--min-mapq N] [--min-baseq N] [--min-expected-frac X] [-f] [--skip-running]
USAGE
}

config=""
sample=""
outdir=""
protocol=""
min_mapq="20"
min_baseq="20"
min_expected_frac="0.8"
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --protocol) protocol="$2"; shift 2 ;;
    --min-mapq) min_mapq="$2"; shift 2 ;;
    --min-baseq) min_baseq="$2"; shift 2 ;;
    --min-expected-frac) min_expected_frac="$2"; shift 2 ;;
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
: "${bamdir:?CONFIG must define bamdir}"

if [ -z "$outdir" ]; then
  outdir="${vcfdir%/}/strand_blacklist"
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
rna_dir_label="${rna_tumor_label:-RNA_TUMOR}"
phased_ext="${strand_blacklist_vcf_extension:-${rna7_phased_vcf_extension:-${phased_vcf_extension:-}}}"
rna_bam_suffix="${strand_blacklist_rna_bam_suffix:-${rna7_smfixed_bam_suffix:-}}"
vep_dir="${strand_blacklist_vep_dir:-${rna_edit_vep_dir:-${mupexi_outdir:-${datadir%/}/mupexi2}}}"
protocol="${protocol:-${strand_blacklist_protocol:-${rna_strand_protocol:-fr-firststrand}}}"

[ -n "$phased_ext" ] || { echo "ERROR: phased VCF extension not found in CONFIG" >&2; exit 1; }
[ -n "$rna_bam_suffix" ] || { echo "ERROR: RNA BAM suffix not found in CONFIG" >&2; exit 1; }

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
inputs=()
patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac
  sid="$(printf '%s\n' "$line" | awk -F'[\t, ]+' '{print $1}')"
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
  bam="${bamdir}/${patient}_${rna_dir_label}/${patient}_${rna_bam_suffix}"
  vep="${vep_dir}/${patient}_vep.vep"
  if [ ! -f "$vcf" ]; then
    echo "[skip] ${patient}: missing phased VCF: $vcf"
    continue
  fi
  if [ ! -f "$bam" ]; then
    echo "[skip] ${patient}: missing RNA BAM: $bam"
    continue
  fi
  if [ ! -f "$vep" ]; then
    echo "[skip] ${patient}: missing VEP file: $vep"
    continue
  fi

  inputs+=("--input" "${patient}=${vcf}=${bam}")
  patients+=("$patient")
done < "$samples"

if [ "${#inputs[@]}" -eq 0 ]; then
  echo "ERROR: no eligible inputs found for strand blacklist" >&2
  exit 1
fi

cohort_support_outfile="${outdir}/cohort.strand_support.tsv"
cohort_blacklist_outfile="${outdir}/cohort.strand_blacklist.tsv"
cohort_uploaded_variation_outfile="${outdir}/cohort.strand_blacklist.uploaded_variation.txt"
tag="cohort"
expected_output="$cohort_blacklist_outfile"
if [ -n "$sample" ]; then
  tag="${patients[0]}"
  expected_output="${outdir}/${tag}.strand_blacklist.tsv"
fi

if [ "$force" != "1" ] && [ -s "$expected_output" ]; then
  echo "[skip] output already exists: ${expected_output} (use -f to overwrite)"
  exit 0
fi

logroot="${outdir}/strand_blacklist.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_strand_blacklist"
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
apply_inputs=()
for item in "${inputs[@]}"; do
  apply_inputs+=("$(printf '%q' "$item")")
done

cohort_args=""
if [ -z "$sample" ]; then
  cohort_args="--cohort-support-outfile $(printf '%q' "$cohort_support_outfile") --cohort-blacklist-outfile $(printf '%q' "$cohort_blacklist_outfile") --cohort-uploaded-variation-outfile $(printf '%q' "$cohort_uploaded_variation_outfile")"
fi

cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${repo_dir}/research/extract_rna_edit_strand_blacklist.py" \\
  ${apply_inputs[*]} \\
  --outdir "$(printf '%q' "$outdir")" \\
  --vep-dir "$(printf '%q' "$vep_dir")" \\
  --protocol "$(printf '%q' "$protocol")" \\
  --min-mapq "$(printf '%q' "$min_mapq")" \\
  --min-baseq "$(printf '%q' "$min_baseq")" \\
  --min-expected-frac "$(printf '%q' "$min_expected_frac")" \\
  ${cohort_args}

if [ ! -s "$(printf '%q' "$expected_output")" ]; then
  echo "ERROR: strand blacklist output missing/empty: $(printf '%q' "$expected_output")" >&2
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
echo "[info] strand blacklist output dir: ${outdir}"
echo ".. logs and reports saved in ${logroot}"
