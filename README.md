# Comprehensive variant calling pipeline for model organisms

This pipeline consists of:

* fastq processing with bbduk.sh
* Read mapping with BWA-mem
* SNP/indel calling with GATK
* Structural variant calling with Meerkat

The pipeline is designed to be applied to model organisms, where reference annotation sets for SNPs/indels may not be available.
