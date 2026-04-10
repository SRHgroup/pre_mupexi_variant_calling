# post_rnadnavar_mupexi_prep (PBS/qsub pipeline)

This repository contains two coordinated modules:

- RNA-editing post-processing: `rna1..rna7`
- Germline calling/filtering: `gdna1..gdna4`

## Layout

- `post_rnadnavar_mupexi_prep/`: RNA step scripts + `run_all.sh`
- `post_rnadnavar_mupexi_prep/rnae_script/`: RNA helper scripts (`rnae1..rnae7`)
- `germline_calling/`: gDNA step scripts + `run_all.sh`
- `germline_calling/scripts/`: helper python
- `examples/CONFIG.example`: example of a CONFIG file, which you write in run_pipeline wrapper, it should have the paths to references, desired file extensions and your data folder structure. 
- `examples/SAMPLES.example`:example of a sample sheet, that you also feed to a run_pipeline wrapper. 
- `bin/check_outputs.sh`: commands that print you a input/output status for a step you want to run.
- `run_all_end_to_end.sh`: run whole rna/gdna end-to-end (works partially at the moment)

## Step names

RNA:
1. `rna1_OnlyRnaVcf.sh` # remove mutations from RNA vcf that exists in DNA-tumor VCF.
2. `rna2_FilterEditSignature.sh` # only keep mutations that follow ADAR and APOBEC3 signature. Adds the specific field in INFO.
3. `rna3_AnnotateKnownSites.sh` # Annotate a specific KNOWN_DB field in the INFO, inditating whether cancer RNA-editing databases have this event. These are be default the only variants later passed to neoantigen preiction.
4. `rna4_SummariseRnaMetrics.sh` # A step needded specifically to default rnadnavar RNA VCF, as it has separate sample tracks for readgroups, which we summarise into one.
5. `rna5_FilterByAfDpAr.sh` # Add optionally QC filters to remove potential artefacts.
6. `rna6_MergeDnaRnaVcfs.sh` # Merge with VCF containing somatic and germline variants.
7. `rna7.0_FixRnaBamReadGroups.sh` (required before `rna7`) # similar to rna4, removes read group tag from reads in RNA BAM, so the phaser can connect reads in BAM with. Make sure you use both BAM files afer GATK preprocessing for both DNA and RNA. 
8. `rna7_GenotypeAndPhaseMergedVcf.sh` # standardises and phases all calls from every sourse of variants into MuPeXi2-ready format. 

gDNA:
1. `gdna1_HaplotypeCaller.sh` # runs germline variant calling
2. `gdna2_FilterGermline.sh` # applies basic QC filter labels
3. `gdna3_SelectVariants.sh` # removes low QC entries 
4. `gdna4_FilterGermlineByAdjacency.sh` # for a given k and somatic/rna-editing VCF preserved only variants adjacent to at least one cancer mutation. Optional but reccomended, as raw germline calls are huge.

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
- `source_rna_mutect2_vcf_extension`, `source_dna_mutect2_vcf_extension` if your upstream rnadnavar VCF names differ from defaults
  - `{patient}` placeholder is supported in these suffixes

Labels:
- default: `DNA_TUMOR`, `RNA_TUMOR`
- UK spelling supported via label vars in `CONFIG`. Make sure you use UK or US spelling everywhere, and save your mental health.

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

Splicing pipeline, step 1: merge STAR junction shards in place:
```bash
make run_splicing_merge_star_sj CONFIG=/path/to/CONFIG
make run_splicing_merge_star_sj CONFIG=/path/to/CONFIG SAMPLE=Pat21
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

If you use the project-folder wrapper (`examples/run_pipeline.sh`), you can check one step per patient or forthe whole cohort:

```bash
./run_pipeline.sh check-step rna5 # runs rna5 script for the whole cohort
./run_pipeline.sh check-step gdna4 SampleX # runs germline variant calling only for SampleX
```

## Notes

- Scripts precheck inputs before `qsub`; broken/missing inputs are reported immediately.
- Duplicate submission guard is enabled per step/sample using stored job IDs + `qstat`.
- `rna7` requires the SM-fixed RNA BAM from `rna7.0` and will fail/skip if that file is missing.
