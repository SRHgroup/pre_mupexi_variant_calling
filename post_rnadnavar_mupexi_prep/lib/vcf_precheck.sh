#!/usr/bin/bash
set -euo pipefail

# Shared precheck helpers for VCF inputs across pipeline steps.
# Default is intentionally tolerant on login-node checks:
# - require file existence
# - require valid gzip if .gz
# - do NOT hard-fail on header grep quirks unless PRECHECK_STRICT_HEADER=1

precheck_vcfgz() {
  local path="$1"
  local tag="${2:-precheck}"
  local hint="${3:-input VCF required}"
  if [ ! -f "$path" ]; then
    echo "[precheck] ${tag}: missing input VCF: $path (${hint})" >&2
    return 1
  fi
  if ! gzip -t "$path" >/dev/null 2>&1; then
    echo "[precheck] ${tag}: invalid gzip VCF: $path" >&2
    return 1
  fi

  if [ "${PRECHECK_STRICT_HEADER:-0}" = "1" ]; then
    if ! zgrep -m1 '^#CHROM' "$path" >/dev/null 2>&1; then
      echo "[precheck] ${tag}: malformed VCF header (missing #CHROM): $path" >&2
      return 1
    fi
  fi
  return 0
}

precheck_vcf_or_vcfgz() {
  local path="$1"
  local tag="${2:-precheck}"
  local hint="${3:-input VCF required}"
  if [ ! -f "$path" ]; then
    echo "[precheck] ${tag}: missing input VCF: $path (${hint})" >&2
    return 1
  fi
  if [[ "$path" = *.gz ]]; then
    precheck_vcfgz "$path" "$tag" "$hint"
    return $?
  fi

  if [ "${PRECHECK_STRICT_HEADER:-0}" = "1" ]; then
    if ! grep -m1 '^#CHROM' "$path" >/dev/null 2>&1; then
      echo "[precheck] ${tag}: malformed VCF header (missing #CHROM): $path" >&2
      return 1
    fi
  fi
  return 0
}

