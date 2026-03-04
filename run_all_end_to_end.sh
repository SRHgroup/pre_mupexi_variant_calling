#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash run_all_end_to_end.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]

Behavior:
- Submits `gdna1..gdna4` as one dependent chain.
- Submits `rna1..rna5` as an independent dependent chain.
- Submits `rna6 -> rna7.0 -> rna7` only after BOTH chains complete successfully.
USAGE
}

force_flag=""
sample_flag=""

while :; do
  case ${1:-} in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c requires path" >&2; exit 1; }
      config=$2
      shift
      ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s requires value" >&2; exit 1; }
      sample_flag="-s $2"
      shift
      ;;
    -f|--force)
      force_flag="-f"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *) break ;;
  esac
  shift
done

[ -n "${config:-}" ] || { usage; exit 1; }
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

submit_chain() {
  local label="$1"
  local initial_depend="$2"
  shift 2
  local steps=("$@")
  local current_depend="$initial_depend"
  local tmp
  tmp="$(mktemp)"

  for step in "${steps[@]}"; do
    echo "[end_to_end] launching ${label}/${step}"
    if [ -n "$current_depend" ]; then
      echo "[end_to_end] dependency: $current_depend"
    fi

    if [[ "$label" == "rna" ]]; then
      # shellcheck disable=SC2086
      if ! QSUB_DEPEND="$current_depend" bash "${repo_root}/post_rnadnavar_mupexi_prep/${step}" -c "$config" $sample_flag $force_flag 2>&1 | tee "$tmp"; then
        rm -f "$tmp"
        return 1
      fi
    else
      # shellcheck disable=SC2086
      if ! QSUB_DEPEND="$current_depend" bash "${repo_root}/germline_calling/${step}" -c "$config" $sample_flag $force_flag 2>&1 | tee "$tmp"; then
        rm -f "$tmp"
        return 1
      fi
    fi

    mapfile -t step_jobids < <(awk -F'jobid=' '/^\[submit\] /{print $2}' "$tmp" | awk '{print $1}' | sed '/^$/d')
    if [ "${#step_jobids[@]}" -gt 0 ]; then
      current_depend="afterok:$(IFS=:; echo "${step_jobids[*]}")"
      echo "[end_to_end] ${label}/${step}: collected ${#step_jobids[@]} job(s)"
    else
      current_depend=""
      echo "[end_to_end] ${label}/${step}: no new jobs submitted"
    fi
  done

  rm -f "$tmp"
  LAST_DEP="$current_depend"
}

extract_ids() {
  local dep="$1"
  dep="${dep#afterok:}"
  printf '%s\n' "$dep"
}

LAST_DEP=""
submit_chain "gdna" "" \
  gdna1_HaplotypeCaller.sh \
  gdna2_FilterGermline.sh \
  gdna3_SelectVariants.sh \
  gdna4_FilterGermlineByAdjacency.sh
dep_gdna="$LAST_DEP"

LAST_DEP=""
submit_chain "rna" "" \
  rna1_OnlyRnaVcf.sh \
  rna2_FilterEditSignature.sh \
  rna3_AnnotateKnownSites.sh \
  rna4_SummariseRnaMetrics.sh \
  rna5_FilterByAfDpAr.sh
dep_rna_pre="$LAST_DEP"

combined_depend=""
if [ -n "$dep_gdna" ] && [ -n "$dep_rna_pre" ]; then
  combined_depend="afterok:$(extract_ids "$dep_gdna"):$(extract_ids "$dep_rna_pre")"
elif [ -n "$dep_gdna" ]; then
  combined_depend="$dep_gdna"
elif [ -n "$dep_rna_pre" ]; then
  combined_depend="$dep_rna_pre"
fi

LAST_DEP=""
submit_chain "rna" "$combined_depend" \
  rna6_MergeDnaRnaVcfs.sh \
  rna7.0_FixRnaBamReadGroups.sh \
  rna7_GenotypeAndPhaseMergedVcf.sh

echo "[end_to_end] submission complete"
