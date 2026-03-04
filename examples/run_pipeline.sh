#!/usr/bin/bash
set -euo pipefail

# Copy this file to your project folder, then edit these two paths:
REPO="/home/projects/SRHgroup/apps/pre_mupexi_variant_calling"
CONFIG="/home/projects/SRHgroup/projects/PemBOv_trial/bin/rnadnavar/CONFIG"

cmd="${1:-help}"
shift || true

force=0
if [[ " ${*:-} " == *" -f "* ]] || [[ " ${*:-} " == *" --force "* ]]; then
  force=1
  set -- "${@/-f/}"
  set -- "${@/--force/}"
fi
force_arg=""
if [ "$force" -eq 1 ]; then
  force_arg="FORCE=1"
fi

run_make() {
  local target="$1"
  local sample="${2:-}"
  if [ -n "$sample" ]; then
    make -C "$REPO" "$target" CONFIG="$CONFIG" SAMPLE="$sample" $force_arg
  else
    make -C "$REPO" "$target" CONFIG="$CONFIG" $force_arg
  fi
}

case "$cmd" in
  rna)       run_make run_all_rna "${1:-}" ;;
  germline)  run_make run_all_germline "${1:-}" ;;
  all)       run_make run_all "${1:-}" ;;
  step)
    step_name="${1:-}"
    sample="${2:-}"
    case "$step_name" in
      rna1|rna2|rna3|rna4|rna5|rna6|rna7|gdna1|gdna2|gdna3|gdna4)
        run_make "run${step_name}" "$sample"
        ;;
      *)
        echo "Unknown step: $step_name" >&2
        exit 1
        ;;
    esac
    ;;
  check)
    sample="${1:-}"
    mode="${2:-all}"
    if [ -n "$sample" ]; then
      make -C "$REPO" check_outputs CONFIG="$CONFIG" SAMPLE="$sample" MODE="$mode"
    else
      make -C "$REPO" check_outputs CONFIG="$CONFIG" MODE="$mode"
    fi
    ;;
  sync)
    git -C "$REPO" pull --ff-only origin main
    ;;
  show-config)
    echo "REPO=$REPO"
    echo "CONFIG=$CONFIG"
    ;;
  *)
    cat <<USAGE
Usage:
  $0 rna [PATIENT] [-f]
  $0 germline [PATIENT] [-f]
  $0 all [PATIENT] [-f]
  $0 step <rna1|rna2|rna3|rna4|rna5|rna6|rna7|gdna1|gdna2|gdna3|gdna4> [PATIENT] [-f]
  $0 check [PATIENT] [all|rna|germline]
  $0 sync
  $0 show-config
USAGE
    ;;
esac
