#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: bash run_all_end_to_end.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [-f]
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

# shellcheck disable=SC2086
bash germline_calling/run_all.sh -c "$config" $sample_flag $force_flag
# shellcheck disable=SC2086
bash post_rnadnavar_mupexi_prep/run_all.sh -c "$config" $sample_flag $force_flag
