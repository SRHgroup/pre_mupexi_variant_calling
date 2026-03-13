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

q_nodes="${strand_blacklist_qsub_nodes:-1}"
q_ppn="${strand_blacklist_qsub_ppn:-2}"
q_mem="${strand_blacklist_qsub_mem:-16gb}"
q_walltime="${strand_blacklist_qsub_walltime:-08:00:00}"
agg_nodes="${strand_blacklist_concat_qsub_nodes:-1}"
agg_ppn="${strand_blacklist_concat_qsub_ppn:-1}"
agg_mem="${strand_blacklist_concat_qsub_mem:-4gb}"
agg_walltime="${strand_blacklist_concat_qsub_walltime:-01:00:00}"

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

declare -a dependency_jobids=()
declare -a aggregate_support_files=()
declare -a aggregate_blacklist_files=()
declare -a aggregate_uploaded_files=()
submitted_count=0
eligible_count=0

logroot="${outdir}/strand_blacklist.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_strand_blacklist"

declare -A seen_patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac
  sid="$(printf '%s\n' "$line" | awk -F'[\t, ]+' '{print $1}')"
  patient="$(sample_base_name "$sid")"
  [ -n "$patient" ] || continue
  [ -n "${seen_patients[$patient]:-}" ] && continue
  seen_patients["$patient"]=1

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

  eligible_count=$((eligible_count + 1))

  support_out="${outdir}/${patient}.strand_support.tsv"
  blacklist_out="${outdir}/${patient}.strand_blacklist.tsv"
  uploaded_out="${outdir}/${patient}.strand_blacklist.uploaded_variation.txt"
  aggregate_support_files+=("$support_out")
  aggregate_blacklist_files+=("$blacklist_out")
  aggregate_uploaded_files+=("$uploaded_out")

  if [ "$force" != "1" ] && [ -s "$support_out" ] && [ -s "$blacklist_out" ]; then
    echo "[skip] ${patient}: outputs already exist (use -f to overwrite)"
    continue
  fi

  job_name="${prefix}.${patient}"
  marker="${logdir}/submitted.${prefix}.${patient}.jobid"
  active_jobid=""
  if [ -f "$marker" ]; then
    prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ]; then
      pbs_state="$(pbs_state_for_jobid "$prev_jobid")"
      if [ "$pbs_state" = "RUNNING" ] || [ "$pbs_state" = "QUEUED" ]; then
        active_jobid="$prev_jobid"
      fi
    fi
  fi
  if [ -z "$active_jobid" ] && command -v qselect >/dev/null 2>&1; then
    active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "${job_name}" 2>/dev/null | head -n1 || true)"
  fi
  if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
    active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="${job_name}" '$4==n {print $1; exit}')"
  fi
  if [ -n "$active_jobid" ]; then
    echo "[skip] ${job_name}: scheduler already has active job ${active_jobid}"
    dependency_jobids+=("$active_jobid")
    continue
  fi

  runscript="${logdir}/run.${patient}.${prefix}.sh"
  cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${repo_dir}/research/extract_rna_edit_strand_blacklist.py" \\
  --input "$(printf '%q' "${patient}=${vcf}=${bam}")" \\
  --outdir "$(printf '%q' "$outdir")" \\
  --vep-dir "$(printf '%q' "$vep_dir")" \\
  --protocol "$(printf '%q' "$protocol")" \\
  --min-mapq "$(printf '%q' "$min_mapq")" \\
  --min-baseq "$(printf '%q' "$min_baseq")" \\
  --min-expected-frac "$(printf '%q' "$min_expected_frac")"

if [ ! -s "$(printf '%q' "$blacklist_out")" ]; then
  echo "ERROR: strand blacklist output missing/empty: $(printf '%q' "$blacklist_out")" >&2
  exit 2
fi
if [ ! -s "$(printf '%q' "$support_out")" ]; then
  echo "ERROR: strand support output missing/empty: $(printf '%q' "$support_out")" >&2
  exit 3
fi
SCRIPT
  chmod +x "$runscript"

  qsub_opts=()
  [ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
  [ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
  qsub_opts+=(-N "${job_name}")
  qsub_opts+=(-o "${repdir}/${prefix}.${patient}.o\$PBS_JOBID")
  qsub_opts+=(-e "${repdir}/${prefix}.${patient}.e\$PBS_JOBID")
  qsub_opts+=(-l "nodes=${q_nodes}:ppn=${q_ppn},mem=${q_mem},walltime=${q_walltime}")

  jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
  printf '%s\n' "$jobid" > "$marker"
  dependency_jobids+=("$jobid")
  submitted_count=$((submitted_count + 1))
  echo "[submit] ${job_name}: jobid=${jobid}"
done < "$samples"

if [ "$eligible_count" -eq 0 ]; then
  echo "ERROR: no eligible inputs found for strand blacklist" >&2
  exit 1
fi

if [ -n "$sample" ]; then
  echo "[info] strand blacklist output dir: ${outdir}"
  echo ".. logs and reports saved in ${logroot}"
  exit 0
fi

cohort_support_outfile="${outdir}/cohort.strand_support.tsv"
cohort_blacklist_outfile="${outdir}/cohort.strand_blacklist.tsv"
cohort_uploaded_variation_outfile="${outdir}/cohort.strand_blacklist.uploaded_variation.txt"
aggregate_job_name="${prefix}.cohort"
aggregate_marker="${logdir}/submitted.${aggregate_job_name}.jobid"

if [ "$force" != "1" ] && [ -s "$cohort_support_outfile" ] && [ -s "$cohort_blacklist_outfile" ] && [ "${#dependency_jobids[@]}" -eq 0 ]; then
  echo "[skip] cohort outputs already exist: ${cohort_blacklist_outfile} (use -f to overwrite)"
  echo "[info] strand blacklist output dir: ${outdir}"
  echo ".. logs and reports saved in ${logroot}"
  exit 0
fi

active_aggregate_jobid=""
if [ -f "$aggregate_marker" ]; then
  prev_jobid="$(head -n1 "$aggregate_marker" 2>/dev/null || true)"
  if [ -n "$prev_jobid" ]; then
    pbs_state="$(pbs_state_for_jobid "$prev_jobid")"
    if [ "$pbs_state" = "RUNNING" ] || [ "$pbs_state" = "QUEUED" ]; then
      active_aggregate_jobid="$prev_jobid"
    fi
  fi
fi
if [ -z "$active_aggregate_jobid" ] && command -v qselect >/dev/null 2>&1; then
  active_aggregate_jobid="$(qselect -u "${USER:-$(whoami)}" -N "${aggregate_job_name}" 2>/dev/null | head -n1 || true)"
fi
if [ -z "$active_aggregate_jobid" ] && command -v qstat >/dev/null 2>&1; then
  active_aggregate_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="${aggregate_job_name}" '$4==n {print $1; exit}')"
fi
if [ -n "$active_aggregate_jobid" ]; then
  echo "[skip] ${aggregate_job_name}: scheduler already has active job ${active_aggregate_jobid}"
  echo "[info] strand blacklist output dir: ${outdir}"
  echo ".. logs and reports saved in ${logroot}"
  exit 0
fi

aggregate_script="${logdir}/run.cohort.${prefix}.sh"
{
  printf '#!/usr/bin/bash\n'
  printf 'set -euo pipefail\n'
  printf 'support_files=('
  for path in "${aggregate_support_files[@]}"; do printf ' %q' "$path"; done
  printf ' )\n'
  printf 'blacklist_files=('
  for path in "${aggregate_blacklist_files[@]}"; do printf ' %q' "$path"; done
  printf ' )\n'
  printf 'uploaded_files=('
  for path in "${aggregate_uploaded_files[@]}"; do printf ' %q' "$path"; done
  printf ' )\n'
  printf 'cohort_support_outfile=%q\n' "$cohort_support_outfile"
  printf 'cohort_blacklist_outfile=%q\n' "$cohort_blacklist_outfile"
  printf 'cohort_uploaded_variation_outfile=%q\n' "$cohort_uploaded_variation_outfile"
  cat <<'SCRIPT'
support_existing=()
for f in "${support_files[@]}"; do
  if [ -s "$f" ]; then
    support_existing+=("$f")
  else
    echo "[warn] missing/empty support file during concat: $f" >&2
  fi
done

blacklist_existing=()
for f in "${blacklist_files[@]}"; do
  if [ -s "$f" ]; then
    blacklist_existing+=("$f")
  else
    echo "[warn] missing/empty blacklist file during concat: $f" >&2
  fi
done

if [ "${#support_existing[@]}" -eq 0 ]; then
  echo "ERROR: no per-patient support files available to concatenate" >&2
  exit 2
fi
if [ "${#blacklist_existing[@]}" -eq 0 ]; then
  echo "ERROR: no per-patient blacklist files available to concatenate" >&2
  exit 3
fi

awk 'FNR==1 && NR!=1 {next} {print}' "${support_existing[@]}" > "$cohort_support_outfile"
awk 'FNR==1 && NR!=1 {next} {print}' "${blacklist_existing[@]}" > "$cohort_blacklist_outfile"

: > "$cohort_uploaded_variation_outfile"
for f in "${uploaded_files[@]}"; do
  [ -f "$f" ] || continue
  cat "$f" >> "$cohort_uploaded_variation_outfile"
done
sort -u -o "$cohort_uploaded_variation_outfile" "$cohort_uploaded_variation_outfile"

if [ ! -s "$cohort_support_outfile" ] || [ ! -s "$cohort_blacklist_outfile" ]; then
  echo "ERROR: cohort strand blacklist concat produced missing/empty outputs" >&2
  exit 4
fi
SCRIPT
} > "$aggregate_script"
chmod +x "$aggregate_script"

qsub_opts=()
[ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
[ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
qsub_opts+=(-N "${aggregate_job_name}")
qsub_opts+=(-o "${repdir}/${aggregate_job_name}.o\$PBS_JOBID")
qsub_opts+=(-e "${repdir}/${aggregate_job_name}.e\$PBS_JOBID")
qsub_opts+=(-l "nodes=${agg_nodes}:ppn=${agg_ppn},mem=${agg_mem},walltime=${agg_walltime}")
if [ "${#dependency_jobids[@]}" -gt 0 ]; then
  dep_string="$(IFS=:; echo "${dependency_jobids[*]}")"
  qsub_opts+=(-W "depend=afterok:${dep_string}")
fi

aggregate_jobid="$(qsub "${qsub_opts[@]}" "$aggregate_script")"
printf '%s\n' "$aggregate_jobid" > "$aggregate_marker"
echo "[submit] ${aggregate_job_name}: jobid=${aggregate_jobid}"
echo "[info] per-patient jobs submitted: ${submitted_count}"
echo "[info] strand blacklist output dir: ${outdir}"
echo ".. logs and reports saved in ${logroot}"
