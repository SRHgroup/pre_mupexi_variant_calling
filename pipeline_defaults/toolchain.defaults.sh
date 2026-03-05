#!/usr/bin/bash
# Default module/toolchain profile for pre_mupexi_variant_calling.
# Override any variable in your project CONFIG or environment.

# Core module versions
module_gatk="gatk/4.5.0.0"
module_bcftools="bcftools/1.23"
module_htslib="htslib/1.23"
module_anaconda="anaconda3/2025.06-1"
module_java="java/17-openjdk"

# Process-specific module sets
modules_rna="ngs tools ${module_htslib} ${module_bcftools} ${module_gatk} ${module_anaconda}"
modules_rna_bamfix="tools ${module_htslib} samtools/1.23"
modules_gdna_hc="tools ngs ${module_java} ${module_gatk}"
modules_gdna_post="tools ngs ${module_anaconda}"
modules_dna_only_merge="ngs tools ${module_htslib} ${module_bcftools} ${module_anaconda}"
modules_dna_only_phase="ngs tools ${module_htslib} ${module_bcftools} ${module_gatk} ${module_anaconda}"
