#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_rna_edit_cluster_jobs.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [--max-distance N] [--min-cluster-size N] [--min-alt-count N] [-f] [--skip-running]
USAGE
}

config=""
sample=""
outdir=""
max_distance=33
min_cluster_size=2
min_alt_count=10
force=0
skip_running=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --max-distance) max_distance="$2"; shift 2 ;;
    --min-cluster-size) min_cluster_size="$2"; shift 2 ;;
    --min-alt-count) min_alt_count="$2"; shift 2 ;;
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

if [ -z "$outdir" ]; then
  outdir="${vcfdir%/}/editing_clusters33"
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"

out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
phased_ext="${rna7_phased_vcf_extension:-${phased_vcf_extension:-rna7.DNAt_DNAn_RNAt_merged_phased.vcf.gz}}"

logroot="${outdir}/rna_edit_cluster_jobs.logs_and_reports"
logdir="${logroot}/logs"
repdir="${logroot}/reports"
mkdir -p "$logdir" "$repdir"

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

  variants_tsv="${outdir}/${patient}.rna_edit_variants.tsv"
  clusters_tsv="${outdir}/${patient}.rna_edit_clusters.tsv"
  plot_prefix="${outdir}/${patient}.rna_edit_cluster_landscape"
  plot_png="${plot_prefix}.png"
  plot_pdf="${plot_prefix}.pdf"

  if [ "$force" != "1" ] && [ -s "$variants_tsv" ] && [ -s "$clusters_tsv" ] && [ -s "$plot_png" ] && [ -s "$plot_pdf" ]; then
    echo "[skip] ${patient}: outputs already exist (use -f to overwrite)"
    continue
  fi

  prefix="research_rna_clusters"
  marker="${logdir}/submitted.${prefix}.${patient}.jobid"
  if [ "$skip_running" = "1" ] && [ -f "$marker" ]; then
    prev_jobid="$(head -n1 "$marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ]; then
      pbs_state="$(pbs_state_for_jobid "$prev_jobid")"
      if [ "$pbs_state" = "RUNNING" ] || [ "$pbs_state" = "QUEUED" ]; then
        echo "[skip-running] ${prefix}.${patient}: active job ${prev_jobid} (${pbs_state})"
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

python3 "${repo_dir}/research/extract_rna_editing_clusters.py" \\
  --input "${patient}=${vcf}" \\
  --out-variants "${variants_tsv}" \\
  --out-clusters "${clusters_tsv}" \\
  --max-distance "${max_distance}" \\
  --min-cluster-size "${min_cluster_size}" \\
  --min-alt-count "${min_alt_count}" \\
  --rna-label "${rna_tumor_label:-RNA_TUMOR}" \\
  --tumor-label "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" \\
  --normal-label "${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"

python3 "${repo_dir}/research/plot_rna_editing_clusters.py" \\
  --variants "${variants_tsv}" \\
  --clusters "${clusters_tsv}" \\
  --out-prefix "${plot_prefix}"
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
