# post_rnadnavar_mupexi_prep (PBS/qsub pipeline)

This repository contains two coordinated modules:

- RNA-editing post-processing: `rna1..rna7`
- Germline calling/filtering: `gdna1..gdna4`

## Layout

- `post_rnadnavar_mupexi_prep/`: RNA step scripts + `run_all.sh`
- `post_rnadnavar_mupexi_prep/rnae_script/`: RNA helper scripts (`rnae1..rnae7`)
- `germline_calling/`: gDNA step scripts + `run_all.sh`
- `germline_calling/scripts/`: helper python
- `examples/CONFIG.example`
- `examples/SAMPLES.example`
- `bin/check_outputs.sh`
- `run_all_end_to_end.sh`
- `Makefile`

## Step names

RNA:
1. `rna1_OnlyRnaVcf.sh`
2. `rna2_FilterEditSignature.sh`
3. `rna3_AnnotateKnownSites.sh`
4. `rna4_SummariseRnaMetrics.sh`
5. `rna5_FilterByAfDpAr.sh`
6. `rna6_MergeDnaRnaVcfs.sh`
7. `rna7.0_FixRnaBamReadGroups.sh` (required before `rna7`)
8. `rna7_GenotypeAndPhaseMergedVcf.sh`

gDNA:
1. `gdna1_HaplotypeCaller.sh`
2. `gdna2_FilterGermline.sh`
3. `gdna3_SelectVariants.sh`
4. `gdna4_FilterGermlineByAdjacency.sh`

## Dependency model

- `rna1..rna5` are independent from `gdna1..gdna4`.
- `rna6`, `rna7.0`, and `rna7` require both branches to be done.
- `run_all_end_to_end.sh` submits exactly that topology:
  - chain A: `gdna1 -> gdna2 -> gdna3 -> gdna4`
  - chain B: `rna1 -> rna2 -> rna3 -> rna4 -> rna5`
  - chain C: `rna6 -> rna7.0 -> rna7`, dependent on A and B completion.

## Requirements

- PBS `qsub` / `qstat`
- `bcftools` / `htslib`
- `gatk`
- `python3`
- `whatshap`
- `samtools`
- `bgzip`

## Configure

```bash
cp examples/CONFIG.example /path/to/project/CONFIG
cp examples/SAMPLES.example /path/to/project/SAMPLES.tsv
```

Edit `CONFIG` and set at minimum:
- `samples`, `vcfdir`, `bamdir`, `knownsites`, `FASTA`, `DICT`, `rnae_scripts`

Labels:
- default: `DNA_TUMOR`, `RNA_TUMOR`
- UK spelling supported via label vars in `CONFIG`

## Run

Single RNA step:
```bash
bash post_rnadnavar_mupexi_prep/rna4_SummariseRnaMetrics.sh -c /path/to/CONFIG -s Pat11
```

Single gDNA step:
```bash
bash germline_calling/gdna2_FilterGermline.sh -c /path/to/CONFIG -s Pat11
```

RNA chain:
```bash
make run_all_rna CONFIG=/path/to/CONFIG
```

gDNA chain:
```bash
make run_all_germline CONFIG=/path/to/CONFIG
```

End-to-end (branching dependencies):
```bash
bash run_all_end_to_end.sh -c /path/to/CONFIG
```

Force recompute for a step:
```bash
make runrna5 CONFIG=/path/to/CONFIG SAMPLE=Pat11 FORCE=1
```

## Check outputs

```bash
bash bin/check_outputs.sh -c /path/to/CONFIG -m all
bash bin/check_outputs.sh -c /path/to/CONFIG -m rna
bash bin/check_outputs.sh -c /path/to/CONFIG -m germline
```

If you use the project-folder wrapper (`examples/run_pipeline.sh`), you can check one step per patient with YES/NO:

```bash
./run_pipeline.sh check-step rna5
./run_pipeline.sh check-step gdna4 01-CH-L
```

## Notes

- Scripts precheck inputs before `qsub`; broken/missing inputs are reported immediately.
- Duplicate submission guard is enabled per step/sample using stored job IDs + `qstat`.
- `rna7` requires the SM-fixed RNA BAM from `rna7.0` and will fail/skip if that file is missing.
- Dataset-specific `CONFIG` files should stay outside git-tracked repo content.
