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
modules_rna="ngs tools ${module_htslib} ${module_bcftools} ${module_java} ${module_gatk} ${module_anaconda}"
modules_rna_bamfix="tools ${module_htslib} samtools/1.23"
modules_gdna_hc="tools ngs ${module_htslib} samtools/1.23 ${module_java} ${module_gatk}"
modules_gdna_post="tools ngs ${module_anaconda}"
modules_dna_only_merge="ngs tools ${module_htslib} ${module_bcftools} ${module_anaconda}"
modules_dna_only_phase="ngs tools ${module_htslib} ${module_bcftools} ${module_java} ${module_gatk} ${module_anaconda}"

# Research/auxiliary jobs
research_python_modules="tools ngs ${module_anaconda}"
splicing_module_prereq="tools"
splicing_python_modules="${module_anaconda}"
splicing_qsub_nodes="1"
splicing_qsub_ppn="2"
splicing_qsub_mem="16gb"
splicing_qsub_walltime="06:00:00"
splicing_spl3_qsub_nodes="${splicing_qsub_nodes}"
splicing_spl3_qsub_ppn="${splicing_qsub_ppn}"
splicing_spl3_qsub_mem="${splicing_qsub_mem}"
splicing_spl3_qsub_walltime="${splicing_qsub_walltime}"
splicing_spl4_qsub_nodes="${splicing_qsub_nodes}"
splicing_spl4_qsub_ppn="${splicing_qsub_ppn}"
splicing_spl4_qsub_mem="24gb"
splicing_spl4_qsub_walltime="08:00:00"
mupexi_modules="tools ngs ${module_anaconda} netmhcpan/4.0a perl/5.36.1 ensembl-tools/90"
