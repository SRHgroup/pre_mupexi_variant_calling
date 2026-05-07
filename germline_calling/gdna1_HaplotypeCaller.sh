#!/usr/bin/bash
set -euo pipefail

# gdna1: call germline variants in DNA normal BAM with GATK HaplotypeCaller
# restricted to somatic mutation loci from DNA tumor and optional RNA tumor VCFs.

usage() {
  cat <<'USAGE'
Usage: bash gdna1_HaplotypeCaller.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
USAGE
}

force=0
while :; do
  case ${1:-} in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c/--config requires a path" >&2; exit 1; }
      config=$2
      shift
      ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s/--sample requires a value" >&2; exit 1; }
      sample=$2
      shift
      ;;
    -f|--force)
      force=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      break
      ;;
  esac
  shift
done

[ -n "${config:-}" ] || { usage; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }
source "$config"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pipeline_defaults="${PIPELINE_DEFAULTS:-${repo_root}/pipeline_defaults/toolchain.defaults.sh}"

: "${samples:?CONFIG must define samples}"
: "${bamdir:?CONFIG must define bamdir}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${FASTA:?CONFIG must define FASTA}"
: "${gdna1_vcf_extension:?CONFIG must define gdna1_vcf_extension}"

dna_normal_label="${dna_normal_label:-DNA_NORMAL}"
dna_bam_suffix="${dna_bam_suffix:-md.bam}"
out_normal_label="${out_dna_normal_label:-${dna_normal_label}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"

sample_base_name() {
  local value="$1"
  local labels=(
    "${dna_normal_label}" "${dna_tumor_label:-DNA_TUMOR}" "${rna_tumor_label:-RNA_TUMOR}"
    "${out_dna_normal_label:-DNA_NORMAL}" "${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}" "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "DNA_NORMAL" "DNA_TUMOR" "DNA_TUMOUR" "RNA_TUMOR" "RNA_TUMOUR" "TUMOR" "TUMOUR"
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

resolve_patient_placeholder() {
  local template="$1"
  local patient="$2"
  printf '%s\n' "${template//\{patient\}/$patient}"
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

detect_source_mutect2_vcfs() {
  local patient="$1"
  local source_rna_ext source_dna_ext
  local rna_vcf_dir dna_vcf_dir srna_pref sdna_pref srna sdna

  source_rna_ext="${source_rna_mutect2_vcf_extension:-${out_rna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
  source_dna_ext="${source_dna_mutect2_vcf_extension:-${out_dna_label}_vs_{patient}_${out_normal_label}.mutect2.filtered.vcf.gz}"
  source_rna_ext="$(resolve_patient_placeholder "$source_rna_ext" "$patient")"
  source_dna_ext="$(resolve_patient_placeholder "$source_dna_ext" "$patient")"
  source_rna_ext="${source_rna_ext%\}}"
  source_dna_ext="${source_dna_ext%\}}"

  rna_vcf_dir="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}"
  dna_vcf_dir="${vcfdir}/${patient}_${out_dna_label}_vs_${patient}_${out_normal_label}"
  srna_pref="${rna_vcf_dir}/${patient}_${source_rna_ext}"
  sdna_pref="${dna_vcf_dir}/${patient}_${source_dna_ext}"

  srna="$(pick_first_existing \
    "$srna_pref" \
    "${rna_vcf_dir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${patient}_${rna_tumor_label:-RNA_TUMOR}_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${patient}_RNA_TUMOR_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${rna_vcf_dir}/${patient}_RNA_TUMOUR_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
  )" || srna=""

  sdna="$(pick_first_existing \
    "$sdna_pref" \
    "${dna_vcf_dir}/${patient}_${out_dna_label}_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${patient}_${dna_tumor_label:-DNA_TUMOR}_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${patient}_DNA_TUMOR_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
    "${dna_vcf_dir}/${patient}_DNA_TUMOUR_vs_${patient}_${out_normal_label}.mutect2.filtered.vcf.gz" \
  )" || sdna=""

  printf '%s\t%s\n' "$sdna" "$srna"
}

prefix=$(basename "${BASH_SOURCE[0]}" .sh)
scriptdir="${vcfdir}/${prefix}.logs_and_reports"
logdir="${scriptdir}/logs"
repdir="${scriptdir}/reports"
mkdir -p "$logdir" "$repdir"
qsub_depend_arg=()
if [ -n "${QSUB_DEPEND:-}" ]; then
  qsub_depend_arg=(-W "depend=${QSUB_DEPEND}")
fi

if [ -z "${sample:-}" ]; then
  echo "Running gdna1 for all patients in $samples"
else
  echo "Running gdna1 only for $sample"
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

  normal_bam="${bamdir}/${name}_${dna_normal_label}/${name}_${dna_normal_label}.${dna_bam_suffix}"
  germline_dir="${vcfdir}/${name}_${out_normal_label}"
  outvcf="${germline_dir}/${name}_${gdna1_vcf_extension}"
  read -r dna_tumor_vcf rna_tumor_vcf <<< "$(detect_source_mutect2_vcfs "$name")"
  sites_bed="${germline_dir}/${name}_gdna1.somatic_sites.bed"
  readlen_tsv="${germline_dir}/${name}_gdna1.read_length.tsv"

  if [ "$force" -eq 0 ] && [ -f "$outvcf" ]; then
    echo "[skip] ${prefix}.${name}: output already exists: $outvcf (use -f to overwrite)"
    continue
  fi
  if [ ! -f "$normal_bam" ]; then
    echo "[precheck] ${prefix}.${name}: missing DNA normal BAM: $normal_bam" >&2
    echo "[skip] ${prefix}.${name}: not submitting qsub due to failed input precheck" >&2
    continue
  fi
  if [ -z "$dna_tumor_vcf" ] && [ -z "$rna_tumor_vcf" ]; then
    echo "[precheck] ${prefix}.${name}: missing both DNA and RNA somatic source VCFs under ${vcfdir}" >&2
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
module load ${modules_gdna_hc:-tools ngs htslib/1.23 samtools/1.23 java/17-openjdk gatk/4.5.0.0}
SCRIPT
    printf 'normal_bam=%q\n' "$normal_bam"
    printf 'germline_dir=%q\n' "$germline_dir"
    printf 'outvcf=%q\n' "$outvcf"
    printf 'FASTA=%q\n' "$FASTA"
    printf 'dna_tumor_vcf=%q\n' "$dna_tumor_vcf"
    printf 'rna_tumor_vcf=%q\n' "$rna_tumor_vcf"
    printf 'sites_bed=%q\n' "$sites_bed"
    printf 'readlen_tsv=%q\n' "$readlen_tsv"
    cat <<'SCRIPT'
if [ ! -f "$normal_bam" ]; then
  echo "ERROR: missing DNA normal BAM: $normal_bam" >&2
  exit 1
fi
if [ -z "${dna_tumor_vcf:-}" ] && [ -z "${rna_tumor_vcf:-}" ]; then
  echo "ERROR: both DNA and RNA somatic source VCFs are missing" >&2
  exit 1
fi

mkdir -p "$germline_dir"

python3 - "$sites_bed" "$dna_tumor_vcf" "$rna_tumor_vcf" <<'PY'
import gzip
import os
import sys

out_path = sys.argv[1]
inputs = [p for p in sys.argv[2:] if p]

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

intervals = set()
for path in inputs:
    if not os.path.exists(path):
        continue
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except Exception:
                continue
            ref = cols[3] if cols[3] not in ("", ".") else "N"
            alts = [a for a in cols[4].split(",") if a and a != "."]
            span = len(ref)
            for alt in alts:
                if alt.startswith("<"):
                    continue
                span = max(span, len(alt))
            if span < 1:
                span = 1
            start = max(0, pos - 1)
            end = start + span
            intervals.add((chrom, start, end))

with open(out_path, "w") as out:
    for chrom, start, end in sorted(intervals, key=lambda x: (x[0], x[1], x[2])):
        out.write(f"{chrom}\t{start}\t{end}\n")
PY

if [ ! -s "$sites_bed" ]; then
  echo "ERROR: somatic sites BED is missing/empty: $sites_bed" >&2
  exit 1
fi

read_length="$(
  samtools view -F 0x900 "$normal_bam" |
  awk '
    NR > 50000 { exit }
    {
      l = length($10)
      if (l > 0) {
        c[l]++
      }
    }
    END {
      best = ""
      bestc = -1
      for (l in c) {
        if (c[l] > bestc || (c[l] == bestc && l + 0 > best + 0)) {
          best = l
          bestc = c[l]
        }
      }
      if (best == "") {
        exit 1
      }
      print best
    }
  '
)"

printf 'metric\tvalue\nread_length\t%s\n' "$read_length" > "$readlen_tsv"
interval_padding="$(awk -F'\t' '$1=="read_length"{print $2; exit}' "$readlen_tsv")"
if ! [[ "$interval_padding" =~ ^[0-9]+$ ]]; then
  echo "ERROR: invalid interval padding extracted from $readlen_tsv: ${interval_padding}" >&2
  exit 1
fi

echo "[info] DNA normal BAM: $normal_bam"
echo "[info] DNA tumor VCF: ${dna_tumor_vcf:-NA}"
echo "[info] RNA tumor VCF: ${rna_tumor_vcf:-NA}"
echo "[info] Somatic sites BED: $sites_bed"
echo "[info] Interval padding (read length): $interval_padding"

gatk HaplotypeCaller \
  -I "$normal_bam" \
  -R "$FASTA" \
  -L "$sites_bed" \
  -ip "$interval_padding" \
  -O "$outvcf" \
  -G StandardAnnotation
SCRIPT
  } > "$runscript"
  chmod +x "$runscript"

  qsub_output="$(qsub -W group_list="${qsub_group:-srhgroup}" -A "${qsub_account:-srhgroup}" -d "$(pwd)" \
    "${qsub_depend_arg[@]}" \
    -l nodes=1:ppn=8,mem=24gb,walltime="00:12:00:00" -r y -N "$job_name" -o "$repdir" -e "$repdir" "$runscript")"
  echo "$qsub_output"
  printf '%s\n' "$qsub_output" > "$submit_marker"
  echo "[submit] ${job_name}: jobid=${qsub_output}"

  echo ".. logs and reports saved in $scriptdir"
  sleep 0.5
done < "$samples"
