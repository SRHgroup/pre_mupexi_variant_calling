# Changelog

All notable changes to this repository should be documented here.

This project follows a simple Keep-a-Changelog style and semantic version tags (`vMAJOR.MINOR.PATCH`).

## [Unreleased]

- Record all pending changes before the next tag.

## [0.1.0] - 2026-03-03

### Added

- Initial import of `post_rnadnavar_mupexi_prep` step scripts (`4.1`-`4.7.1`) and helper scripts.
- Repository structure with `post_rnadnavar_mupexi_prep/`, `post_rnadnavar_mupexi_prep/rnae_script/`, `examples/`, `bin/`.
- `run_all.sh` wrapper, `Makefile` step targets, and `bin/check_outputs.sh` validation utility.
- Documentation for setup, execution, labels, outputs, and HPC clone/update workflow.

### Changed

- Standardized bash shebangs to `#!/usr/bin/bash`.
- Enforced `set -euo pipefail` in step scripts.
- Hardened runscript generation with quoted heredocs + explicit variable assignment.
- Improved sample suffix detection flexibility for `TUMOR`/`TUMOUR` naming.

### Notes for future updates

- Add a new dated section under `Unreleased` first.
- When releasing, move entries into a new version section and create a git tag.
- Keep behavior-impacting changes explicit (inputs/outputs, defaults, thresholds, labels).
