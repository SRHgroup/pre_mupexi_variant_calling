#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash splicing/run_merge_star_sj_shards.sh -c CONFIG [-s SAMPLE_OR_PATIENT] [--root STAR_DIR] [-f] [--dry-run]

Behavior:
- Recursively finds STAR shard files such as *.0001.SJ.out.tab under the STAR root
- Merges shard groups into one STAR-style output next to the shards
- Writes outputs named like Pat21_RNA_TUMOUR.SJ.out.tab
USAGE
}

config=""
sample=""
root_override=""
force=0
dry_run=0

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
    -f|--force)
      force=1
      shift
      ;;
    --dry-run)
      dry_run=1
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
script_path="${repo_root}/splicing/merge_star_sj_shards.py"

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

cmd=(python3 "$script_path" --root "$star_root")
if [ -n "$sample" ]; then
  cmd+=(--sample-filter "$sample")
fi
if [ "$force" -eq 1 ]; then
  cmd+=(--force)
fi
if [ "$dry_run" -eq 1 ]; then
  cmd+=(--dry-run)
fi

printf '[info] STAR root: %s\n' "$star_root"
"${cmd[@]}"
