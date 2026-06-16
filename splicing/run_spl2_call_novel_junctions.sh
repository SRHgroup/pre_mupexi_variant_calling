#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash splicing/run_spl2_call_novel_junctions.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [--root STAR_DIR] [--gtf GTF] [--min-unique-reads N] [-f] [--dry-run] [--include-noncanonical] [--keep-non-protein-coding]

Behavior:
- Builds known splice junctions from CONFIG:GTF or --gtf
- Scans merged STAR files such as *.SJ.out.tab under the STAR root
- Writes *.spl2.novel_junctions.tsv next to each merged STAR file
- Keeps SSNIP's tumor support threshold by default: unique junction reads >= 10
USAGE
}

config=""
sample=""
root_override=""
gtf_override=""
min_unique_reads="10"
force=0
dry_run=0
include_noncanonical=0
keep_non_protein_coding=0

while [ $# -gt 0 ]; do
  case "${1:-}" in
    -c|--config)
      config="${2:-}"
      shift 2
      ;;
    -s|--sample)
      sample="${2:-}"
      shift 2
      ;;
    --root)
      root_override="${2:-}"
      shift 2
      ;;
    --gtf)
      gtf_override="${2:-}"
      shift 2
      ;;
    --min-unique-reads)
      min_unique_reads="${2:-}"
      shift 2
      ;;
    -f|--force)
      force=1
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    --include-noncanonical)
      include_noncanonical=1
      shift
      ;;
    --keep-non-protein-coding)
      keep_non_protein_coding=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[ -n "$config" ] || { usage >&2; exit 1; }
[ -f "$config" ] || { echo "ERROR: config not found: $config" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
script_path="${repo_root}/splicing/call_novel_junctions.py"

if [ -n "$root_override" ]; then
  star_root="$root_override"
elif [ -n "${stardir:-}" ]; then
  star_root="$stardir"
elif [ -n "${bamdir:-}" ]; then
  preprocessing_root="$(dirname "$bamdir")"
  star_root="${preprocessing_root}/star"
else
  echo "ERROR: CONFIG must define stardir or bamdir, or pass --root" >&2
  exit 1
fi

if [ -n "$gtf_override" ]; then
  gtf_path="$gtf_override"
elif [ -n "${GTF:-}" ]; then
  gtf_path="$GTF"
else
  echo "ERROR: CONFIG must define GTF, or pass --gtf" >&2
  exit 1
fi

[ -f "$gtf_path" ] || { echo "ERROR: missing GTF: $gtf_path" >&2; exit 1; }

cmd=(python3 "$script_path" --star-root "$star_root" --gtf "$gtf_path" --min-unique-reads "$min_unique_reads")
if [ -n "$sample" ]; then
  cmd+=(--sample-filter "$sample")
fi
if [ "$force" -eq 1 ]; then
  cmd+=(--force)
fi
if [ "$dry_run" -eq 1 ]; then
  cmd+=(--dry-run)
fi
if [ "$include_noncanonical" -eq 1 ]; then
  cmd+=(--include-noncanonical)
fi
if [ "$keep_non_protein_coding" -eq 1 ]; then
  cmd+=(--keep-non-protein-coding)
fi

printf '[spl2] STAR root: %s\n' "$star_root"
printf '[spl2] GTF: %s\n' "$gtf_path"
"${cmd[@]}"
