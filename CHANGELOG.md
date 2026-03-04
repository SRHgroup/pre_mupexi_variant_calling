# Changelog

All notable changes to this repository should be documented here.

This project follows a simple Keep-a-Changelog style and semantic version tags (`vMAJOR.MINOR.PATCH`).

## [Unreleased]

### Added

- New `germline_calling/` module with:
  - `gdna1_HaplotypeCaller.sh`
  - `gdna2_FilterGermline.sh`
  - `gdna3_SelectVariants.sh`
  - `gdna4_FilterGermlineByAdjacency.sh`
  - `scripts/filter_germline_by_proximity.py`
  - `run_all.sh`
- End-to-end wrapper `run_all_end_to_end.sh`.
- Legacy germline reformat scripts kept in `germline_calling/legacy/`.

### Changed

- Renamed step scripts to semantic names:
  - RNA: `rna1`..`rna7`
  - gDNA: `gdna1`..`gdna4`
- Renamed CONFIG suffix variables to match step names (for example `rna5_qced_vcf_extension`, `gdna4_vcf_extension`).
- Updated `Makefile` targets to `runrna*` / `rungdna*` (legacy numeric target aliases retained).
- Updated end-to-end submission flow so `rna1..rna5` and `gdna1..gdna4` submit independently, while `rna6..rna7` wait for both branches.
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
