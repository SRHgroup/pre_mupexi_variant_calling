# Backlog

## Deferred RNA-Editing QC Improvements

### 1. Add proper strandedness control as `rna5.1`
- Insert a new RNA-edit QC step between `rna5` and `rna6`.
- Goal: strand-aware filtering before RNA edits are merged into the phased DNA/RNA VCF.
- Intended behavior:
  - use transcript strand plus RNA BAM read orientation
  - retain only transcript-consistent RNA-edit support for stranded libraries
  - emit a blacklist/filter output and feed the filtered RNA-edit set into `rna6`
- Motivation:
  - current post hoc blacklist generation is useful for already processed cohorts, but the correct long-term implementation belongs in the core pipeline before merge/phasing

### 2. Add APOBEC3 motif-based confidence labeling
- Replace the current "known DB" style confidence concept for APOBEC3 with motif-based labeling.
- Motivation:
  - there are effectively no trusted known-database matches for APOBEC3-style RNA editing
  - motif context is more informative for confidence than database overlap in this case
- Reference to use when implementing:
  - https://www.sciencedirect.com/science/article/pii/S0022283624004844

### 3. Add homopolymer artifact control
- Add a QC filter/annotation for variants falling in or near homopolymer runs such as `AAAAAA` or `CCCCCC`.
- Motivation:
  - these regions can generate sequencing/alignment artifacts that look like RNA-edit candidates
- Intended behavior:
  - annotate homopolymer context
  - optionally downgrade confidence or filter such sites in the RNA-edit workflow
