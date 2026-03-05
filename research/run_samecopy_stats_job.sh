#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_samecopy_stats_job.sh -c CONFIG [-o OUTFILE] [--window N] [-s PATIENT] [-f] [--skip-running]
USAGE
}

config=""
outfile=""
window=33
sample=""
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -o|--outfile) outfile="$2"; shift 2 ;;
    --window) window="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
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

if [ -z "$outfile" ]; then
  stats_dir="${vcfdir}/samecopy_stats"
  mkdir -p "$stats_dir"
  if [ -n "$sample" ]; then
    outfile="${stats_dir}/${sample}.samecopy_stats_with_rna.tsv"
  else
    outfile="${stats_dir}/cohort.samecopy_stats_with_rna.tsv"
  fi
fi

mkdir -p "$(dirname "$outfile")"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
research_python_modules="${research_python_modules:-tools ngs anaconda3/2025.06-1}"

out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"

germline_ext="${samecopy_germline_vcf_extension:-${ndna_vcf:-${gdna4_vcf_extension:-${gdna3_vcf_extension:-${output_extension_30:-${output_extension_202:-2.0.2.Filtered.Selected.HaplotypeCaller.vcf}}}}}}"
src_dna_ext="${source_dna_mutect2_vcf_extension:-${out_dna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
src_rna_ext="${source_rna_mutect2_vcf_extension:-${out_rna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
merged_ext="${rna7_phased_vcf_extension:-${phased_vcf_extension:-4.7_DNAt_DNAn_RNAt_merged_phased.vcf.gz}}"
tumor_col_label="${out_dna_label:-TUMOR}"

sanitize_mutect2_suffix() {
  local s="$1"
  # remove accidental trailing template braces from bad configs
  s="${s%\}}"
  # if duplicated tail is present, keep only first full extension occurrence
  if [[ "$s" == *".mutect2.filtered.vcf.gz"* ]]; then
    s="${s%%.mutect2.filtered.vcf.gz*}.mutect2.filtered.vcf.gz"
  elif [[ "$s" == *".mutect2.filtered.vcf"* ]]; then
    s="${s%%.mutect2.filtered.vcf*}.mutect2.filtered.vcf"
  fi
  printf '%s\n' "$s"
}

src_dna_ext="$(sanitize_mutect2_suffix "$src_dna_ext")"
src_rna_ext="$(sanitize_mutect2_suffix "$src_rna_ext")"

logroot="${vcfdir}/research_samecopy_stats_with_rna.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

prefix="research_samecopy_stats"
tag="cohort"
[ -n "$sample" ] && tag="$sample"
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

if [ -f "$marker" ]; then
  prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
  if [ -n "$prev_jobid" ]; then
    st="$(pbs_state_for_jobid "$prev_jobid")"
    if [ "$st" = "RUNNING" ] || [ "$st" = "QUEUED" ]; then
      echo "[skip] ${prefix}.${tag}: active job ${prev_jobid} (${st})"
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

if [ "$force" != "1" ] && [ -s "$outfile" ]; then
  echo "[skip] output already exists: $outfile (use -f to overwrite)"
  exit 0
fi

runscript="${logdir}/run.${tag}.${prefix}.sh"
cat > "$runscript" <<SCRIPT
#!/usr/bin/bash
set -euo pipefail
if [ -n "\${PIPELINE_DEFAULTS:-}" ] && [ -f "\$PIPELINE_DEFAULTS" ]; then
  # shellcheck disable=SC1090
  source "\$PIPELINE_DEFAULTS"
fi
module load ${research_python_modules}

python3 "${repo_dir}/research/vcf_samecopy_stats_with_rna.py" \\
  --samples "${samples}" \\
  --vcfdir "${vcfdir}" \\
  --out-normal-label "${out_normal_label}" \\
  --out-dna-label "${out_dna_label}" \\
  --out-rna-label "${out_rna_label}" \\
  --gdna3-ext "${germline_ext}" \\
  --source-dna-ext "${src_dna_ext}" \\
  --source-rna-ext "${src_rna_ext}" \\
  --merged-ext "${merged_ext}" \\
  --tumor-col-label "${tumor_col_label}" \\
  --window "${window}" \\
  --outfile "${outfile}" \\
  ${sample:+--patient "${sample}"}
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
echo ".. logs and reports saved in ${logroot}"
