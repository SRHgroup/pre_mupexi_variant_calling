#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/run_mupexi_jobs.sh -c CONFIG [-s PATIENT] [-o OUTDIR] [--run-fusions] [--fusion-only] [--run-splicing] [--splicing-only] [--hla HLA_STRING] [--expr EXPR_TSV] [--fusion FUSION_ARRIBA_TSV] [--splicing SPL4_TSV] [--nodes N] [--ppn N] [--mem SIZE] [--walltime HH:MM:SS] [-f] [--skip-running]
USAGE
}

config=""
sample=""
outdir=""
run_fusions=0
fusion_only=0
run_splicing=0
splicing_only=0
force=0
skip_running=0
cli_hla=""
cli_expr=""
cli_fusion=""
cli_splicing=""
cli_nodes=""
cli_ppn=""
cli_mem=""
cli_walltime=""

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config) config="$2"; shift 2 ;;
    -s|--sample) sample="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --run-fusions) run_fusions=1; shift ;;
    --fusion-only) fusion_only=1; run_fusions=1; shift ;;
    --run-splicing) run_splicing=1; shift ;;
    --splicing-only) splicing_only=1; run_splicing=1; shift ;;
    --hla) cli_hla="$2"; shift 2 ;;
    --expr) cli_expr="$2"; shift 2 ;;
    --fusion) cli_fusion="$2"; shift 2 ;;
    --splicing) cli_splicing="$2"; shift 2 ;;
    --nodes) cli_nodes="$2"; shift 2 ;;
    --ppn) cli_ppn="$2"; shift 2 ;;
    --mem) cli_mem="$2"; shift 2 ;;
    --walltime) cli_walltime="$2"; shift 2 ;;
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

# Compatibility/fallback config keys
mupexi2_repo="${mupexi2_repo:-${mupexi_repo:-${mupexi_repo_dir:-/home/projects/SRHgroup/apps/mupexi2}}}"
mupexi_netmhc_config="${mupexi_netmhc_config:-${netmhc_config:-${mupexi_netmhcpan_config:-}}}"

: "${samples:?CONFIG must define samples}"
: "${vcfdir:?CONFIG must define vcfdir}"
: "${mupexi2_repo:?CONFIG must define mupexi2_repo}"
: "${mupexi_netmhc_config:?CONFIG must define mupexi_netmhc_config}"
: "${hladir:?CONFIG must define hladir}"

if [ -z "$outdir" ]; then
  outdir="${mupexi_outdir:-${vcfdir}/mupexi2}"
fi
mkdir -p "$outdir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

out_normal_label="${out_dna_normal_label:-${dna_normal_label:-DNA_NORMAL}}"
out_dna_label="${out_dna_tumor_label:-${dna_tumor_label:-DNA_TUMOR}}"
out_rna_label="${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
rna_tumor_sample_name="${mupexi_tumor_sample:-${rna7_signal_sample_label:-TUMOR}}"
dna_tumor_sample_name="${mupexi_dna_only_tumor_sample:-${rna_tumor_sample_name}}"
normal_sample_name="${mupexi_normal_sample:-DNA_NORMAL}"
peptide_lengths="${mupexi_peptide_lengths:-9-11}"
parallel_k="${mupexi_parallel_k:-true}"
enable_germlines="${mupexi_enable_germlines:-true}"
enable_superpeptides="${mupexi_enable_superpeptides:-true}"
enable_rna_edit="${mupexi_enable_rna_edit:-true}"
phased_ext="${rna7_phased_vcf_extension:-${phased_vcf_extension:-}}"
dna_only_phased_ext="${dna_only_phased_vcf_extension:-}"
hld_ext="${mupexi_hla_extension:-${output_extension_14:-1.4.RunStatBootstrapMean.Rstat.txt}}"
hld_direct_ext="${mupexi_hla_direct_extension:-hla1.tab}"
expr_dir="${kaldir:-}"
expr_ext="${output_extension_14:-1.4.RunStatBootstrapMean.Rstat.txt}"
fus_dir="${fus_dir:-${fusdir:-}}"
splicing_root="${mupexi_splicing_outdir:-${splicing_outdir:-${datadir:-}/splicing}}"
q_nodes="${cli_nodes:-${mupexi_qsub_nodes:-1}}"
q_mem="${cli_mem:-${mupexi_qsub_mem:-24gb}}"
q_walltime="${cli_walltime:-${mupexi_qsub_walltime:-24:00:00}}"
junction_only=0
if [ "$fusion_only" = "1" ] || [ "$splicing_only" = "1" ]; then
  junction_only=1
fi
if [ "$junction_only" != "1" ]; then
  if [ -z "$phased_ext" ] && [ -z "$dna_only_phased_ext" ]; then
    echo "ERROR: missing phased VCF extension in CONFIG (rna7_phased_vcf_extension/phased_vcf_extension or dna_only_phased_vcf_extension)" >&2
    exit 1
  fi
fi
if [ "$fusion_only" = "1" ] || [ "$splicing_only" = "1" ]; then
  # Force SNV-derived mutation sources off for junction-only runs.
  enable_germlines="false"
  enable_superpeptides="false"
  enable_rna_edit="false"
fi

count_k_values() {
  local spec="$1"
  awk -v s="$spec" '
    BEGIN{
      n=split(s,a,",")
      for(i=1;i<=n;i++){
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",a[i])
        if(a[i]=="") continue
        if(a[i] ~ /^[0-9]+-[0-9]+$/){
          split(a[i],r,"-")
          lo=r[1]+0; hi=r[2]+0
          if(lo>hi){ t=lo; lo=hi; hi=t }
          for(k=lo;k<=hi;k++) seen[k]=1
        } else if(a[i] ~ /^[0-9]+$/){
          seen[a[i]+0]=1
        }
      }
      c=0
      for(k in seen) c++
      if(c<1) c=1
      print c
    }'
}

k_count="$(count_k_values "$peptide_lengths")"
default_ppn="$k_count"
if [[ "${parallel_k,,}" == "false" || "${parallel_k,,}" == "no" || "${parallel_k}" == "0" ]]; then
  default_ppn="1"
fi
q_ppn="${cli_ppn:-${mupexi_qsub_ppn:-$default_ppn}}"

resolve_patient_placeholder() {
  local template="$1"
  local patient="$2"
  local rna_label="${rna_tumor_label:-RNA_TUMOR}"
  local out_rna_label_value="${out_rna_tumor_label:-$rna_label}"
  local tumor_tag_value="${tumor_tag:-TUMOR}"
  local sample_id="${patient}_${out_rna_label_value}"
  local out="$template"
  out="${out//\{patient\}/$patient}"
  out="${out//\{sample\}/$sample_id}"
  out="${out//\{rna_tumor_label\}/$rna_label}"
  out="${out//\{out_rna_tumor_label\}/$out_rna_label_value}"
  out="${out//\{tumor_tag\}/$tumor_tag_value}"
  printf '%s\n' "$out"
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

find_normal_sample_name() {
  local patient="$1"
  local out=""
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    case "$line" in [[:space:]]*'#'*) continue ;; esac
    sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
    base="$(sample_base_name "$sid")"
    if [ "$base" != "$patient" ]; then
      continue
    fi
    if [[ "$sid" == *"_${dna_normal_label:-DNA_NORMAL}" ]] || [[ "$sid" == *"_${out_dna_normal_label:-DNA_NORMAL}" ]] || [[ "$sid" == *"_DNA_NORMAL" ]] || [[ "$sid" == *"_N" ]]; then
      out="$sid"
      break
    fi
  done < "$samples"
  printf '%s\n' "$out"
}

patient_has_rna_sample() {
  local patient="$1"
  local sid base label
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    case "$line" in [[:space:]]*'#'*) continue ;; esac
    sid="$(printf '%s\n' "$line" | awk -F'[,\t ]+' '{print $1}')"
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

extract_hla_from_file() {
  local path="$1"
  [ -f "$path" ] || return 1
  awk '
    BEGIN { OFS=""; n=0 }
    {
      line=$0
      while (match(line, /HLA-[A-Za-z0-9]+\*[0-9]+:[0-9]+/)) {
        h=substr(line, RSTART, RLENGTH)
        gsub(/\*/, "", h)
        if (!(h in seen)) { seen[h]=1; arr[++n]=h }
        line=substr(line, RSTART+RLENGTH)
      }
      while (match(line, /HLA-[A-Za-z0-9]+:[0-9]+:[0-9]+/)) {
        h=substr(line, RSTART, RLENGTH)
        if (!(h in seen)) { seen[h]=1; arr[++n]=h }
        line=substr(line, RSTART+RLENGTH)
      }
      while (match(line, /HLA-[A-Za-z0-9]+:[0-9]+/)) {
        h=substr(line, RSTART, RLENGTH)
        if (!(h in seen)) { seen[h]=1; arr[++n]=h }
        line=substr(line, RSTART+RLENGTH)
      }
    }
    END {
      for (i=1;i<=n;i++) {
        if (i>1) printf ","
        printf "%s", arr[i]
      }
      printf "\n"
    }
  ' "$path"
}

extract_hla_from_optitype_file() {
  local path="$1"
  [ -f "$path" ] || return 1
  awk -F'\t' '
    function trim(s) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
      return s
    }
    function add_hla(raw, h) {
      raw = trim(raw)
      if (raw == "" || raw == "." || raw == "NA") return
      h = raw
      if (h !~ /^HLA-/) h = "HLA-" h
      gsub(/\*/, "", h)
      if (!(h in seen)) { seen[h]=1; arr[++n]=h }
    }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        key = trim($i)
        if (key != "") hdr[key] = i
      }
      next
    }
    NR >= 2 {
      add_hla($(hdr["A1"]))
      add_hla($(hdr["A2"]))
      add_hla($(hdr["B1"]))
      add_hla($(hdr["B2"]))
      add_hla($(hdr["C1"]))
      add_hla($(hdr["C2"]))
      exit
    }
    END {
      for (i = 1; i <= n; i++) {
        if (i > 1) printf ","
        printf "%s", arr[i]
      }
      printf "\n"
    }
  ' "$path"
}

find_nfcore_optitype_hla_file() {
  local patient="$1"
  local dirs=(
    "${hladir}/optitype/${patient}_NORMAL"
    "${hladir}/optitype/${patient}_DNA_NORMAL"
    "${hladir}/optitype/${patient}_N"
    "${hladir}/nf-core-hlatyping/optitype/${patient}_NORMAL"
    "${hladir}/nf-core-hlatyping/optitype/${patient}_DNA_NORMAL"
    "${hladir}/nf-core-hlatyping/optitype/${patient}_N"
  )
  local dir candidate
  for dir in "${dirs[@]}"; do
    [ -d "$dir" ] || continue
    for candidate in \
      "$dir"/*result*.tsv \
      "$dir"/*Result*.tsv \
      "$dir"/*.tsv \
      "$dir"/*.txt \
      "$dir"/*.csv; do
      [ -f "$candidate" ] || continue
      printf '%s\n' "$candidate"
      return 0
    done
  done
  return 1
}

find_expression_file() {
  local patient="$1"
  local candidates=()
  [ -n "${expr_dir:-}" ] || return 1
  [ -n "${expr_ext:-}" ] || return 1
  candidates+=(
    "${expr_dir}/${patient}_${expr_ext}"
    "${expr_dir}/${patient}_${mupexi_tumor_sample:-${rna7_signal_sample_label:-TUMOR}}_${expr_ext}"
    "${expr_dir}/${patient}_${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}_${expr_ext}"
    "${expr_dir}/${patient}_${rna_tumor_label:-RNA_TUMOR}_${expr_ext}"
    "${expr_dir}/${patient}_TUMOR_${expr_ext}"
    "${expr_dir}/${patient}_RNA_TUMOR_${expr_ext}"
    "${expr_dir}/${patient}_RNA_TUMOUR_${expr_ext}"
  )
  local p
  for p in "${candidates[@]}"; do
    if [ -f "$p" ]; then
      printf '%s\n' "$p"
      return 0
    fi
  done
  return 1
}

find_splicing_file() {
  local patient="$1"
  local root="${splicing_root:-}"
  [ -n "$root" ] || return 1
  [ -d "$root" ] || return 1

  local labels=(
    "${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"
    "${rna_tumor_label:-RNA_TUMOR}"
    "RNA_${tumor_tag:-TUMOR}"
    "RNA_TUMOR"
    "RNA_TUMOUR"
  )
  local seen_labels="" label sample_id candidate
  for label in "${labels[@]}"; do
    [ -n "$label" ] || continue
    if printf '%s\n' "$seen_labels" | grep -Fxq "$label"; then
      continue
    fi
    seen_labels="${seen_labels}
${label}"
    sample_id="${patient}_${label}"
    for candidate in \
      "${root}/${sample_id}/${sample_id}.spl4.neojunctions.tsv" \
      "${root}/${sample_id}/${sample_id}_spl4.neojunctions.tsv" \
      "${root}/${sample_id}.spl4.neojunctions.tsv" \
      "${root}/${sample_id}_spl4.neojunctions.tsv" \
      "${root}/${patient}/${sample_id}.spl4.neojunctions.tsv" \
      "${root}/${patient}/${sample_id}_spl4.neojunctions.tsv"; do
      if [ -f "$candidate" ]; then
        printf '%s\n' "$candidate"
        return 0
      fi
    done
  done

  find "$root" -maxdepth 3 -type f \( \
    -name "${patient}*.spl4.neojunctions.tsv" -o \
    -name "${patient}*_spl4.neojunctions.tsv" \
  \) | sort | head -n 1
}

patient_mupexi_output_exists() {
  local patient="$1"
  if [ "$junction_only" = "1" ]; then
    if [ "$run_splicing" = "1" ] && [ "$run_fusions" = "1" ]; then
      [ -e "${outdir}/${patient}_neojunctions.mupexi" ]
      return $?
    fi
    if [ "$run_splicing" = "1" ]; then
      [ -e "${outdir}/${patient}_neospl.mupexi" ]
      return $?
    fi
    if [ "$run_fusions" = "1" ]; then
      [ -e "${outdir}/${patient}_fus.mupexi" ]
      return $?
    fi
  fi
  find "$outdir" -maxdepth 1 \( -type f -o -type d \) -name "${patient}*" | grep -q .
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

  patient_mode="dna_only"
  patient_enable_rna_edit="false"
  patient_tumor_sample_name="${dna_tumor_sample_name}"
  patient_expr_allowed=0
  patient_fusions_allowed=0
  patient_splicing_allowed=0
  if patient_has_rna_sample "$patient"; then
    patient_mode="dna_rna"
    patient_enable_rna_edit="${enable_rna_edit}"
    patient_tumor_sample_name="${rna_tumor_sample_name}"
    patient_expr_allowed=1
    patient_fusions_allowed=1
    patient_splicing_allowed=1
  fi

  vcf=""
  if [ "$junction_only" != "1" ]; then
    if [ "$patient_mode" = "dna_rna" ]; then
      if [ -z "$phased_ext" ]; then
        echo "[skip] ${patient}: RNA present in SAMPLES but CONFIG lacks rna7/phased VCF extension"
        continue
      fi
      vcf="${vcfdir}/${patient}_${out_rna_label}_vs_${patient}_${out_normal_label}/${patient}_${phased_ext}"
    else
      if [ -z "$dna_only_phased_ext" ]; then
        echo "[skip] ${patient}: DNA-only patient but CONFIG lacks dna_only_phased_vcf_extension"
        continue
      fi
      vcf="${vcfdir}/${patient}_${out_dna_label}_vs_${patient}_${out_normal_label}/${patient}_${dna_only_phased_ext}"
    fi
    if [ ! -f "$vcf" ]; then
      echo "[skip] ${patient}: missing ${patient_mode} phased VCF: $vcf"
      continue
    fi
  fi

  expr=""
  if [ "$patient_expr_allowed" = "1" ]; then
    if [ -n "$sample" ] && [ -n "$cli_expr" ]; then
      expr="$cli_expr"
    elif [ -n "$expr_dir" ]; then
      expr="$(find_expression_file "$patient" || true)"
    fi
    if [ -n "$expr" ] && [ ! -f "$expr" ]; then
      echo "[warn] ${patient}: expression file not found, running MuPeXI without -e: $expr"
      expr=""
    elif [ -z "$expr" ]; then
      echo "[warn] ${patient}: no expression configured/found, running MuPeXI without -e"
    fi
  elif [ -n "$sample" ] && [ -n "$cli_expr" ]; then
    echo "[warn] ${patient}: CLI expression override ignored for DNA-only patient"
  fi

  hla=""
  if [ -n "$sample" ] && [ -n "$cli_hla" ]; then
    hla="$cli_hla"
  else
    hla_file=""
    hla_file_type=""
    # Primary expected layout: Pat101_hla1.tab
    direct_hla="${hladir}/${patient}_${hld_direct_ext}"
    if [ -f "$direct_hla" ]; then
      hla_file="$direct_hla"
      hla_file_type="legacy"
    fi
    normal_name="$(find_normal_sample_name "$patient")"
    if [ -n "$normal_name" ]; then
      cand1="${hladir}/${normal_name}_${hld_ext}"
      if [ -z "$hla_file" ] && [ -f "$cand1" ]; then
        hla_file="$cand1"
        hla_file_type="legacy"
      fi
    fi
    # Legacy fallback logic:
    #   name -> name without last char + "1_N", and special H1 -> H2_N case.
    if [ -z "$hla_file" ]; then
      base="${patient%?}"
      cand2="${hladir}/${base}1_N_${hld_ext}"
      if [ -f "$cand2" ]; then
        hla_file="$cand2"
        hla_file_type="legacy"
      fi
    fi
    if [ -z "$hla_file" ] && [[ "$patient" =~ ^(H1)$ ]]; then
      base="${patient%?}"
      cand3="${hladir}/${base}2_N_${hld_ext}"
      if [ -f "$cand3" ]; then
        hla_file="$cand3"
        hla_file_type="legacy"
      fi
    fi
    if [ -z "$hla_file" ]; then
      if nfcore_hla="$(find_nfcore_optitype_hla_file "$patient" || true)"; then
        if [ -n "$nfcore_hla" ] && [ -f "$nfcore_hla" ]; then
          hla_file="$nfcore_hla"
          hla_file_type="optitype"
        fi
      fi
    fi
    # Optional map fallback if explicitly configured.
    if [ -z "$hla_file" ] && [ -n "${mupexi_hla_map_tsv:-}" ]; then
      hla="$(lookup_map_value "$mupexi_hla_map_tsv" "$patient" || true)"
    elif [ -n "$hla_file" ]; then
      if [ "$hla_file_type" = "optitype" ]; then
        hla="$(extract_hla_from_optitype_file "$hla_file" || true)"
      else
        hla="$(extract_hla_from_file "$hla_file" || true)"
      fi
    fi
  fi
  if [ -z "$hla" ]; then
    echo "[skip] ${patient}: missing HLA (expected ${hladir}/${patient}_${hld_direct_ext}; optional override --hla)"
    continue
  fi

  fusion_path=""
  if [ "$run_fusions" = "1" ]; then
    if [ "$patient_fusions_allowed" != "1" ]; then
      if [ "$fusion_only" = "1" ]; then
        echo "[skip] ${patient}: --fusion-only requested for DNA-only patient"
        continue
      else
        echo "[warn] ${patient}: --run-fusions requested for DNA-only patient; running SNVs only"
      fi
    else
      if [ -n "$sample" ] && [ -n "$cli_fusion" ]; then
        fusion_path="$cli_fusion"
      elif [ -n "${mupexi_fusion_arriba_map_tsv:-}" ]; then
        fusion_path="$(lookup_map_value "$mupexi_fusion_arriba_map_tsv" "$patient" || true)"
      elif [ -n "${mupexi_fusion_arriba_template:-}" ]; then
        fusion_path="$(resolve_patient_placeholder "$mupexi_fusion_arriba_template" "$patient")"
      elif [ -n "${fus_dir:-}" ]; then
        for p in \
          "${fus_dir}/${patient}.fusion_arriba.tsv" \
          "${fus_dir}/${patient}.fusions_arriba.tsv" \
          "${fus_dir}/${patient}_fusion_arriba.tsv" \
          "${fus_dir}/${patient}_fusions_arriba.tsv" \
          "${fus_dir}/${patient}/fusion_arriba.tsv" \
          "${fus_dir}/${patient}/fusions_arriba.tsv" \
          "${fus_dir}/${patient}/${patient}.fusion_arriba.tsv" \
          "${fus_dir}/${patient}/${patient}.fusions_arriba.tsv"; do
          if [ -f "$p" ]; then
            fusion_path="$p"
            break
          fi
        done
      fi
      if [ -n "$fusion_path" ] && [ ! -f "$fusion_path" ]; then
        if [ "$fusion_only" = "1" ]; then
          echo "[skip] ${patient}: --fusion-only requested but fusion_arriba file missing"
          continue
        else
          echo "[warn] ${patient}: --run-fusions requested but fusion_arriba file missing; running without -z"
          fusion_path=""
        fi
      fi
      if [ "$fusion_only" = "1" ] && [ -z "$fusion_path" ]; then
        echo "[skip] ${patient}: --fusion-only requested but no fusion file resolved"
        continue
      fi
    fi
  fi

  splicing_path=""
  if [ "$run_splicing" = "1" ]; then
    if [ "$patient_splicing_allowed" != "1" ]; then
      if [ "$splicing_only" = "1" ]; then
        echo "[skip] ${patient}: --splicing-only requested for DNA-only patient"
        continue
      else
        echo "[warn] ${patient}: --run-splicing requested for DNA-only patient; running without -S"
      fi
    else
      if [ -n "$sample" ] && [ -n "$cli_splicing" ]; then
        splicing_path="$cli_splicing"
      elif [ -n "${mupexi_splicing_map_tsv:-}" ]; then
        splicing_path="$(lookup_map_value "$mupexi_splicing_map_tsv" "$patient" || true)"
      elif [ -n "${mupexi_splicing_template:-}" ]; then
        splicing_path="$(resolve_patient_placeholder "$mupexi_splicing_template" "$patient")"
      else
        splicing_path="$(find_splicing_file "$patient" || true)"
      fi
      if [ -n "$splicing_path" ] && [ ! -f "$splicing_path" ]; then
        if [ "$splicing_only" = "1" ]; then
          echo "[skip] ${patient}: --splicing-only requested but spl4 neojunction file missing: $splicing_path"
          continue
        else
          echo "[warn] ${patient}: --run-splicing requested but spl4 neojunction file missing; running without -S: $splicing_path"
          splicing_path=""
        fi
      fi
      if [ "$splicing_only" = "1" ] && [ -z "$splicing_path" ]; then
        echo "[skip] ${patient}: --splicing-only requested but no spl4 neojunction file resolved"
        continue
      fi
    fi
  fi

  # Heuristic skip if patient-prefixed outputs already exist.
  if [ "$force" != "1" ] && patient_mupexi_output_exists "$patient"; then
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

cmd=(
  python3 -m mupexi2.cli
  --germlines "${enable_germlines}"
  --superpeptides "${enable_superpeptides}"
  --rna-edit "${patient_enable_rna_edit}"
  --parallel-k "${parallel_k}"
  -l "${peptide_lengths}"
  -a "${hla}"
  -t -f -n
  -c "${mupexi_netmhc_config}"
  -p "${patient}"
  -d "${outdir}"
)

if [ "$junction_only" != "1" ]; then
  cmd+=(
    -v "${vcf}"
    --vcf-type merged
    --tumor-sample "${patient_tumor_sample_name}"
    --normal-sample "${normal_sample_name}"
  )
fi

if [ -n "${expr}" ]; then
  cmd+=(-e "${expr}")
fi
if [ -n "${fusion_path}" ]; then
  cmd+=(-z "${fusion_path}")
fi
if [ -n "${splicing_path}" ]; then
  cmd+=(-S "${splicing_path}")
fi

"\${cmd[@]}"
SCRIPT
  chmod +x "$runscript"

  qsub_opts=()
  [ -n "${qsub_group:-}" ] && qsub_opts+=(-W "group_list=${qsub_group}")
  [ -n "${qsub_account:-}" ] && qsub_opts+=(-A "${qsub_account}")
  qsub_opts+=(-l "nodes=${q_nodes}:ppn=${q_ppn},mem=${q_mem},walltime=${q_walltime}")
  qsub_opts+=(-N "${prefix}.${patient}")
  qsub_opts+=(-o "${repdir}/${prefix}.${patient}.o\$PBS_JOBID")
  qsub_opts+=(-e "${repdir}/${prefix}.${patient}.e\$PBS_JOBID")

  jobid="$(qsub "${qsub_opts[@]}" "$runscript")"
  printf '%s\n' "$jobid" > "$marker"
  echo "[submit] ${prefix}.${patient}: mode=${patient_mode} jobid=${jobid}"
done < "$samples"

echo ".. logs and reports saved in ${logroot}"
