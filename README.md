# post_rnadnavar_mupexi_prep (PBS/qsub RNA-DNA extension pipeline)

This repository packages the `post_rnadnavar_mupexi_prep/` bash pipeline (steps `4.1` to `4.7.1`) as a versioned, reusable Git project.

## Repository layout

- `post_rnadnavar_mupexi_prep/`: step scripts and `run_all.sh`
- `post_rnadnavar_mupexi_prep/rnae_script/`: helper scripts called by steps
- `germline_calling/`: germline calling/filtering scripts and helper python
- `examples/CONFIG.example`: template config
- `examples/SAMPLES.example`: template sample list
- `bin/check_outputs.sh`: expected-output checker
- `Makefile`: convenience targets

## Pipeline steps

1. `4.1_OnlyRnaVcf.sh`: remove RNA variants where CHROM+POS exists in DNA tumour VCF.
2. `4.2_FilterEditSignature.sh`: split multi-allelic records and keep ADAR/APOBEC3 signatures.
3. `4.3_AnnotateKnownSites.sh`: annotate INFO field with known RNA-editing DB hits.
4. `4.4_SummariseRnaMetrics.sh`: collapse multiple RNA tumour samples into one standardized RNA sample.
5. `4.5_FilterByAfDpAr.sh`: apply DP/AF/ALT-read thresholds.
6. `4.6_MergeDnaRnaVcfs.sh`: merge DNA normal germline + DNA tumour somatic + RNA editing VCFs.
7. `4.7.0_FixRnaBamReadGroups.sh` (optional): rewrite RNA BAM SM tags.
8. `4.7.1_GenotypeAndPhaseMergedVcf.sh`: genotype union alleles and phase (whatshap).

## Germline steps

1. `2.0_HaplotypeCaller.sh`: call germline variants from DNA normal BAM.
2. `2.0.1_FilterGermline.sh`: QC filter (DP/QD and clustered sites).
3. `2.0.2_SelectVariants.sh`: keep non-filtered calls.
4. `3.0_FilterGermlineByAdjacency.sh`: optional proximity filter to DNA somatic sites (and optionally RNA sites).

## Requirements

Cluster environment with PBS `qsub` and typical modules/tools:

- `bcftools` / `htslib`
- `gatk` (HaplotypeCaller)
- `python3`
- `whatshap`
- `samtools`
- `bgzip`

Scripts load modules inside each generated runscript; adjust `module load ...` lines for your HPC stack.

## Configure

1. Create your project config:

```bash
cp examples/CONFIG.example /path/to/your/project/CONFIG
cp examples/SAMPLES.example /path/to/your/project/SAMPLES
```

2. Edit `CONFIG` and set at minimum:

- `samples`
- `vcfdir`
- `bamdir`
- `knownsites`
- `FASTA`
- `DICT`
- `rnae_scripts` (point to this repo clone: `<clone>/post_rnadnavar_mupexi_prep/rnae_script`)
- `ndna_vcf` (germline suffix used by merge step 4.6)

3. Label controls:

- Default is US spelling (`DNA_TUMOR`, `RNA_TUMOR`).
- To switch to UK spelling, set:

```bash
dna_tumor_label="DNA_TUMOUR"
rna_tumor_label="RNA_TUMOUR"
out_dna_tumor_label="DNA_TUMOUR"
out_rna_tumor_label="RNA_TUMOUR"
```

Detection is suffix-based (not hardcoded exact full sample IDs), including optional RNA suffixes like `.0001`.

4. `SAMPLES` structure:

- One row per biological sample (for example DNA normal, DNA tumor, RNA tumor).
- First column must include patient ID + sample-type suffix, e.g. `Pat1_DNA_NORMAL`, `Pat1_DNA_TUMOR`, `Pat1_RNA_TUMOR`.
- The pipeline derives unique patient IDs by stripping configured labels, so each patient is processed once even with 3 rows.
- Columns after the first are currently ignored by the step scripts (you can still keep FASTQ paths and sample_type columns for compatibility with other pipelines).

## Run

From repository `post_rnadnavar_mupexi_prep/` directory:

```bash
bash 4.1_OnlyRnaVcf.sh -c /path/to/CONFIG
bash 4.2_FilterEditSignature.sh -c /path/to/CONFIG
bash 4.3_AnnotateKnownSites.sh -c /path/to/CONFIG
bash 4.4_SummariseRnaMetrics.sh -c /path/to/CONFIG
bash 4.5_FilterByAfDpAr.sh -c /path/to/CONFIG
bash 4.6_MergeDnaRnaVcfs.sh -c /path/to/CONFIG
bash 4.7.1_GenotypeAndPhaseMergedVcf.sh -c /path/to/CONFIG
```

Single-sample run:

```bash
bash 4.4_SummariseRnaMetrics.sh -c /path/to/CONFIG -s Pat11
```

`-s` accepts either full sample ID (for example `Pat11_RNA_TUMOR`) or patient ID (`Pat11`).

Force recompute:

```bash
bash 4.4_SummariseRnaMetrics.sh -c /path/to/CONFIG -s Pat11 -f
```

Convenience commands:

```bash
make run4.1 CONFIG=/path/to/CONFIG
make run_all_rna CONFIG=/path/to/CONFIG
make run_all_germline CONFIG=/path/to/CONFIG
make run_all CONFIG=/path/to/CONFIG
make check_outputs CONFIG=/path/to/CONFIG MODE=all
make check_outputs CONFIG=/path/to/CONFIG MODE=rna
make check_outputs CONFIG=/path/to/CONFIG MODE=germline
```

Full end-to-end wrapper:

```bash
bash run_all_end_to_end.sh -c /path/to/CONFIG
```

## Output and logs

- Per-patient output pattern:
  - Germline outputs: `${vcfdir}/${patient}_DNA_NORMAL/`
  - `${vcfdir}/${patient}_RNA_TUMOR_vs_${patient}_DNA_NORMAL/`
  - (label strings reflect your CONFIG output labels)
- PBS runscripts/log dirs per step:
  - `${vcfdir}/4.X_*.logs_and_reports/logs`
  - `${vcfdir}/4.X_*.logs_and_reports/reports`

## Validation quick checks

List missing expected outputs:

```bash
bash bin/check_outputs.sh -c /path/to/CONFIG
```

Count phased outputs produced:

```bash
find "$(awk -F'=' '/^vcfdir=/{gsub(/"/,"",$2); print $2}' /path/to/CONFIG)" -name '*.4.7_DNAt_DNAn_RNAt_merged_phased.vcf.gz' | wc -l
```

## HPC clone/update workflow (run from other folders)

Recommended pattern:

1. Clone once in a stable tools path (example):

```bash
git clone <repo-url> /hpc/tools/post_rnadnavar_mupexi_prep
```

2. Keep dataset-specific `CONFIG` files outside the repo (one per dataset).
3. Run any step from any working directory by passing absolute `-c /dataset/.../CONFIG`.
4. Update pipeline safely:

```bash
cd /hpc/tools/post_rnadnavar_mupexi_prep
git fetch --tags
git checkout v0.1.0    # or a newer tag/branch
```

5. Test updates on one dataset/single sample before full reruns (use `-s` and optionally `-f`).

This gives reproducible versions (`git checkout <tag>`) while allowing per-dataset parameter changes in external CONFIG files.

## Notes

- All numbered scripts in both pipelines use the same `CONFIG` and `SAMPLES`.
- `SAMPLES` can contain multiple rows per patient (DNA normal, DNA tumor, RNA tumor); scripts process unique patient IDs by suffix stripping.
- Legacy reformat scripts from the older germline workflow are preserved under `germline_calling/legacy/` but are not required for the current whatshap-based path.
