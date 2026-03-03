# Changelog

All notable changes to this repository should be documented here.

This project follows a simple Keep-a-Changelog style and semantic version tags (`vMAJOR.MINOR.PATCH`).

## [Unreleased]

### Added

- New `germline_calling/` module with:
  - `2.0_HaplotypeCaller.sh`
  - `2.0.1_FilterGermline.sh`
  - `2.0.2_SelectVariants.sh`
  - `3.0_FilterGermlineByAdjacency.sh`
  - `scripts/filter_germline_by_proximity.py`
  - `run_all.sh`
- End-to-end wrapper `run_all_end_to_end.sh`.
- Legacy germline reformat scripts kept in `germline_calling/legacy/`.

### Changed

- Unified both pipelines on patient-centric `SAMPLES` parsing with configurable labels (`DNA_NORMAL`, `DNA_TUMOR`, `RNA_TUMOR` by default).
- Updated examples and docs for shared CONFIG/SAMPLES and full germline+RNA workflow.
- Extended `Makefile` with germline and end-to-end targets.

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
