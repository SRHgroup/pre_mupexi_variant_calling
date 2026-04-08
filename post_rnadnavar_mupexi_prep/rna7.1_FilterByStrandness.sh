#!/usr/bin/bash
set -euo pipefail

# rna7.1: post hoc VEP-only strandedness filter on the final phased VCF.
# This is a convenience recovery step. Preferred method is rna5.1 before merge/phasing.

usage() {
  cat <<'USAGE'
Usage: bash rna7.1_FilterByStrandness.sh -c CONFIG [-s SAMPLE] [-f]

Post hoc VEP-only filter on the final phased VCF.
Preferred method remains rna5.1_FilterByStrandness.sh before merge/phasing.
USAGE
}

force=0
while :; do
  case ${1:-} in
    -c|--config)
      if [ -n "${2:-}" ]; then config=$2; shift; else echo "ERROR: -c/--config requires a path" >&2; exit 1; fi ;;
    -s|--sample)
      if [ -n "${2:-}" ]; then sample=$2; shift; else echo "ERROR: -s/--sample requires a sample name" >&2; exit 1; fi ;;
    -f|--force) force=1 ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
  shift
done

[ -n "${config:-}" ] || { echo "ERROR: Config file needed. Use -c </path/to/CONFIG>." >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: Cannot find config file: $config" >&2; exit 1; }
source "$config"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pipeline_defaults="${PIPELINE_DEFAULTS:-${repo_root}/pipeline_defaults/toolchain.defaults.sh}"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${rnae_scripts:?CONFIG must define rnae_scripts}"
: "${rna7_phased_vcf_extension:?CONFIG must define rna7_phased_vcf_extension}"

default_posthoc_vep_dir="${mupexi_outdir:-}"
if [ -z "$default_posthoc_vep_dir" ] && [ -n "${datadir:-}" ]; then
  default_posthoc_vep_dir="${datadir%/}/mupexi2"
fi
vep_dir="${rna7_posthoc_vep_dir:-${rna_edit_vep_dir:-$default_posthoc_vep_dir}}"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label:-DNA_NORMAL}" "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}"
    "${out_dna_normal_label:-DNA_NORMAL}" "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
  )
  local label
  for label in "${labels[@]}"; do value="${value%_${label}}"; done
  printf '%s\n' "$value"
}

sample_is_requested() {
  local sample_id="$1"
  local patient_id="$2"
  local requested="${sample:-}"
  [ -z "$requested" ] || [ "$sample_id" = "$requested" ] || [ "$patient_id" = "$requested" ]
}

pick_first_existing() {
  local path
  for path in "$@"; do
    [ -n "$path" ] || continue
    if [ -f "$path" ]; then
      printf '%s\n' "$path"
      return 0
    fi
  done
  return 1
}

# shellcheck disable=SC1090
source "${repo_root}/post_rnadnavar_mupexi_prep/lib/vcf_precheck.sh"

if [ -z "${sample:-}" ]; then
  echo "Running rna7.1 for all samples in $samples"
else
  echo "Running rna7.1 only for $sample"
fi
echo "[warn] rna7.1 is a post hoc VEP-only convenience filter. Preferred method is rna5.1_FilterByStrandness.sh."

prefix=$(basename "${BASH_SOURCE[0]}" .sh)
scriptdir="${vcfdir}/${prefix}.logs_and_reports"
logdir="${scriptdir}/logs"
repdir="${scriptdir}/reports"
mkdir -p "$logdir" "$repdir"
qsub_depend_arg=()
if [ -n "${QSUB_DEPEND:-}" ]; then
  qsub_depend_arg=(-W "depend=${QSUB_DEPEND}")
fi

declare -A seen_patients=()
while IFS= read -r line; do
  [ -n "$line" ] || continue
  case "$line" in [[:space:]]*'#'*) continue ;; esac

  sample_name=$(printf '%s\n' "$line" | awk -F'[,	 ]+' '{print $1}')
  name=$(sample_base_name "$sample_name")
  [ -n "$name" ] || continue
  sample_is_requested "$sample_name" "$name" || continue
  [[ -n "${seen_patients[$name]:-}" ]] && continue
  seen_patients["$name"]=1

  out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
  out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
  outdir_only="${vcfdir}/${name}_${out_rna_label}_vs_${name}_${out_normal_label}"

  phased_vcf="${outdir_only}/${name}_${rna7_phased_vcf_extension}"
  backup_vcf="${phased_vcf%.vcf.gz}.pre_rna7.1.vcf.gz"
  keep_tsv="${phased_vcf%.vcf.gz}.rna7.1.strand_keep.tsv"
  blacklist_tsv="${phased_vcf%.vcf.gz}.rna7.1.strand_blacklist.tsv"
  stats_file="${phased_vcf%.vcf.gz}.rna7.1.strand_filter_stats.tsv"
  vep_file="$(pick_first_existing "${vep_dir}/${name}_vep.vep" "${vep_dir}/${name}_vep.vep.gz")" || vep_file=""
  source_vcf="${backup_vcf}"
  if [ ! -f "$source_vcf" ]; then
    source_vcf="${phased_vcf}"
  fi

  if [ "$force" -eq 0 ] && [ -f "$stats_file" ] && [ -s "$stats_file" ] && [ -f "$phased_vcf" ] && [ -s "$phased_vcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $stats_file (use -f to overwrite)"
    continue
  fi
  if ! precheck_vcfgz "$source_vcf" "${prefix}.${name}" "run rna7 first"; then
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi
  if [ -z "$vep_file" ]; then
    echo "[precheck] ${prefix}.${name}: missing VEP file in ${vep_dir} (${name}_vep.vep[.gz])" >&2
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi

  job_name="${prefix}.${name}"
  if [ "${SKIP_RUNNING:-0}" = "1" ]; then
    active_jobid=""
    if command -v qselect >/dev/null 2>&1; then
      active_jobid="$(qselect -u "${USER:-$(whoami)}" -N "$job_name" 2>/dev/null | head -n1 || true)"
    fi
    if [ -z "$active_jobid" ] && command -v qstat >/dev/null 2>&1; then
      active_jobid="$(qstat -u "${USER:-$(whoami)}" 2>/dev/null | awk -v n="$job_name" '$4==n {print $1; exit}')"
    fi
    if [ -n "$active_jobid" ]; then
      echo "[skip] ${job_name}: job already active in scheduler: ${active_jobid}"
      continue
    fi
  fi
  submit_marker="${logdir}/submitted.${job_name}.jobid"
  if [ -f "$submit_marker" ]; then
    prev_jobid="$(head -n1 "$submit_marker" 2>/dev/null || true)"
    if [ -n "$prev_jobid" ] && command -v qstat >/dev/null 2>&1 && qstat "$prev_jobid" >/dev/null 2>&1; then
      echo "[skip] ${job_name}: job already queued/running: ${prev_jobid}"
      continue
    fi
  fi

  runscript="${logdir}/run.${name}.${prefix}.sh"
  {
    printf 'export PIPELINE_DEFAULTS=%q\n' "$pipeline_defaults"
    cat <<'SCRIPT'
#!/usr/bin/bash
set -euo pipefail
if [ -n "${PIPELINE_DEFAULTS:-}" ] && [ -f "$PIPELINE_DEFAULTS" ]; then
  source "$PIPELINE_DEFAULTS"
fi
module load ${modules_rna:-ngs tools htslib/1.23 bcftools/1.23 anaconda3/2025.06-1}
SCRIPT
    printf 'source_vcf=%q\n' "$source_vcf"
    printf 'phased_vcf=%q\n' "$phased_vcf"
    printf 'backup_vcf=%q\n' "$backup_vcf"
    printf 'keep_tsv=%q\n' "$keep_tsv"
    printf 'blacklist_tsv=%q\n' "$blacklist_tsv"
    printf 'stats_file=%q\n' "$stats_file"
    printf 'vep_file=%q\n' "$vep_file"
    printf 'patient_name=%q\n' "$name"
    printf 'rnae_scripts=%q\n' "$rnae_scripts"
    cat <<'SCRIPT'
echo "[warn] rna7.1 is a post hoc VEP-only convenience filter. Preferred method is rna5.1_FilterByStrandness.sh."

if [ ! -f "$source_vcf" ]; then
  echo "ERROR: missing phased VCF (run rna7 first): $source_vcf" >&2
  exit 1
fi
if [ ! -f "$vep_file" ]; then
  echo "ERROR: missing VEP file: $vep_file" >&2
  exit 1
fi

if [ ! -f "$backup_vcf" ]; then
  cp "$source_vcf" "$backup_vcf"
  bcftools index -f -t "$backup_vcf"
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

split_vcf="$tmpdir/input.split.vcf.gz"
filtered_raw="$tmpdir/output.filtered.vcf"
filtered_vcfgz="$tmpdir/output.filtered.vcf.gz"

bcftools norm -m -any -Oz -o "$split_vcf" "$backup_vcf"
bcftools index -f -t "$split_vcf"

python3 "$rnae_scripts/rnae7_1_filter_by_vep_strandedness.py" \
  -i "$split_vcf" \
  -o "$filtered_raw" \
  --vep "$vep_file" \
  --patient "$patient_name" \
  --keep-tsv "$keep_tsv" \
  --blacklist-tsv "$blacklist_tsv" \
  --stats "$stats_file"

bgzip -f -c "$filtered_raw" > "$filtered_vcfgz"

if ! bcftools view -h "$filtered_vcfgz" | grep -q '^#CHROM'; then
  echo "ERROR: rna7.1 produced malformed VCF (missing #CHROM): $filtered_vcfgz" >&2
  exit 1
fi

mv "$filtered_vcfgz" "$phased_vcf"
bcftools index -f -t "$phased_vcf"
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=2,mem=16gb,walltime="00:06:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
