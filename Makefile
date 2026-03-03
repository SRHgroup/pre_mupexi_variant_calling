SHELL := /usr/bin/bash

CONFIG ?=
SAMPLE ?=
FORCE ?=

SAMPLE_FLAG := $(if $(SAMPLE),-s $(SAMPLE),)
FORCE_FLAG := $(if $(filter 1 true yes,$(FORCE)),-f,)

check-config:
	@if [[ -z "$(CONFIG)" ]]; then echo "Set CONFIG=/path/to/CONFIG"; exit 1; fi

run4.1: check-config
	cd rnadnavar && bash 4.1_OnlyRnaVcf.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run4.2: check-config
	cd rnadnavar && bash 4.2_FilterEditSignature.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run4.3: check-config
	cd rnadnavar && bash 4.3_AnnotateKnownSites.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run4.4: check-config
	cd rnadnavar && bash 4.4_SummariseRnaMetrics.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run4.5: check-config
	cd rnadnavar && bash 4.5_FilterByAfDpAr.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run4.6: check-config
	cd rnadnavar && bash 4.6_MergeDnaRnaVcfs.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run4.7: check-config
	cd rnadnavar && bash 4.7.1_GenotypeAndPhaseMergedVcf.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

run_all: check-config
	cd rnadnavar && bash run_all.sh -c "$(CONFIG)" $(SAMPLE_FLAG) $(FORCE_FLAG)

check_outputs: check-config
	bash bin/check_outputs.sh -c "$(CONFIG)" $(SAMPLE_FLAG)
