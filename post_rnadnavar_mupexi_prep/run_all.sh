#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash run_all.sh -c CONFIG [-s SAMPLE] [-f]
USAGE
}

force_flag=""
sample_flag=""

while :; do
  case ${1:-} in
    -c|--config)
      [ -n "${2:-}" ] || { echo "ERROR: -c requires a path" >&2; exit 1; }
      config=$2
      shift
      ;;
    -s|--sample)
      [ -n "${2:-}" ] || { echo "ERROR: -s requires a sample" >&2; exit 1; }
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
    *)
      break
      ;;
  esac
  shift
done

[ -n "${config:-}" ] || { usage; exit 1; }

current_depend="${QSUB_DEPEND:-}"

run_step_collect() {
  local step="$1"
  local tmp
  tmp="$(mktemp)"

  echo "[run_all] launching $step"
  if [ -n "$current_depend" ]; then
    echo "[run_all] dependency: $current_depend"
  fi

  # shellcheck disable=SC2086
  if ! QSUB_DEPEND="$current_depend" bash "$step" -c "$config" $sample_flag $force_flag 2>&1 | tee "$tmp"; then
    rm -f "$tmp"
    return 1
  fi

  mapfile -t step_jobids < <(awk -F'jobid=' '/^\[submit\] /{print $2}' "$tmp" | awk '{print $1}' | sed '/^$/d')
  rm -f "$tmp"

  if [ "${#step_jobids[@]}" -gt 0 ]; then
    current_depend="afterok:$(IFS=:; echo "${step_jobids[*]}")"
    echo "[run_all] collected ${#step_jobids[@]} job(s) from $step"
  else
    current_depend=""
    echo "[run_all] no new jobs submitted by $step"
  fi
}

for step in \
  rna1_OnlyRnaVcf.sh \
  rna2_FilterEditSignature.sh \
  rna3_AnnotateKnownSites.sh \
  rna4_SummariseRnaMetrics.sh \
  rna5_FilterByAfDpAr.sh \
  rna6_MergeDnaRnaVcfs.sh \
  rna7.0_FixRnaBamReadGroups.sh \
  rna7_GenotypeAndPhaseMergedVcf.sh \
  rna8_ReshapeVcf.sh
  do
  run_step_collect "$step"
  done
