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

for step in \
  4.1_OnlyRnaVcf.sh \
  4.2_FilterEditSignature.sh \
  4.3_AnnotateKnownSites.sh \
  4.4_SummariseRnaMetrics.sh \
  4.5_FilterByAfDpAr.sh \
  4.6_MergeDnaRnaVcfs.sh \
  4.7.1_GenotypeAndPhaseMergedVcf.sh
  do
  echo "[run_all] launching $step"
  # shellcheck disable=SC2086
  bash "$step" -c "$config" $sample_flag $force_flag
  done
