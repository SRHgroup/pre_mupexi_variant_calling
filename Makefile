SHELL := /usr/bin/bash

CONFIG ?=
SAMPLE ?=
FORCE ?=
MODE ?= all
OUTDIR ?=
MAX_DISTANCE ?= 33
MIN_CLUSTER_SIZE ?= 2
MIN_ALT_COUNT ?= 10
WINDOW ?= 33
OUTFILE ?=

SAMPLE_FLAG := $(if $(SAMPLE),-s $(SAMPLE),)
FORCE_FLAG := $(if $(filter 1 true yes,$(FORCE)),-f,)
SKIP_RUNNING_FLAG := $(if $(filter 1 true yes,$(SKIP_RUNNING)),--skip-running,)
OUTDIR_FLAG := $(if $(OUTDIR),-o $(OUTDIR),)
OUTFILE_FLAG := $(if $(OUTFILE),-o $(OUTFILE),)

check-config:
	@if [[ -z "$(CONFIG)" ]]; then echo "Set CONFIG=/path/to/CONFIG"; exit 1; fi

runrna1: check-config
	cd post_rnadnavar_mupexi_prep && bash rna1_OnlyRnaVcf.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna2: check-config
	cd post_rnadnavar_mupexi_prep && bash rna2_FilterEditSignature.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna3: check-config
	cd post_rnadnavar_mupexi_prep && bash rna3_AnnotateKnownSites.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna4: check-config
	cd post_rnadnavar_mupexi_prep && bash rna4_SummariseRnaMetrics.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna5: check-config
	cd post_rnadnavar_mupexi_prep && bash rna5_FilterByAfDpAr.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna6: check-config
	cd post_rnadnavar_mupexi_prep && bash rna6_MergeDnaRnaVcfs.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna7.0: check-config
	cd post_rnadnavar_mupexi_prep && bash rna7.0_FixRnaBamReadGroups.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

runrna7: check-config
	cd post_rnadnavar_mupexi_prep && bash rna7_GenotypeAndPhaseMergedVcf.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

rungdna1: check-config
	cd germline_calling && bash gdna1_HaplotypeCaller.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

rungdna2: check-config
	cd germline_calling && bash gdna2_FilterGermline.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

rungdna3: check-config
	cd germline_calling && bash gdna3_SelectVariants.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

rungdna4: check-config
	cd germline_calling && bash gdna4_FilterGermlineByAdjacency.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run_all_rna: check-config
	cd post_rnadnavar_mupexi_prep && bash run_all.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run_all_germline: check-config
	cd germline_calling && bash run_all.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run_all: check-config
	bash run_all_end_to_end.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run_dna_only: check-config
	cd post_rnadnavar_mupexi_prep && bash dna_only_merge_and_phase.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run_research_rna_clusters: check-config
	cd research && bash run_rna_edit_cluster_jobs.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(OUTDIR_FLAG) --max-distance "$(MAX_DISTANCE)" --min-cluster-size "$(MIN_CLUSTER_SIZE)" --min-alt-count "$(MIN_ALT_COUNT)" $(FORCE_FLAG) $(SKIP_RUNNING_FLAG)

run_research_samecopy_stats: check-config
	cd research && bash run_samecopy_stats_job.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(OUTFILE_FLAG) --window "$(WINDOW)" $(FORCE_FLAG) $(SKIP_RUNNING_FLAG)

check_outputs: check-config
	bash bin/check_outputs.sh -c "$(CONFIG)" $(SAMPLE_FLAG) -m "$(MODE)"

# Backward-compatible target aliases
run4.1: runrna1
run4.2: runrna2
run4.3: runrna3
run4.4: runrna4
run4.5: runrna5
run4.6: runrna6
run4.7.0: runrna7.0
run4.7: runrna7
run2.0: rungdna1
run2.0.1: rungdna2
run2.0.2: rungdna3
run3.0: rungdna4
