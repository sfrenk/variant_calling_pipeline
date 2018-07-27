# Comprehensive variant calling pipeline for model organisms

This pipeline consists of:

* fastq processing with bbduk.sh
* Read mapping with BWA-mem
* SNP/indel calling with GATK
* Functional annotation of SNPS/indels with snpEff
* Structural variant calling with SvABA and Meerkat
* Copy number estimation of repeat elements
* CNV calling with CNVnator
* Identification of novel transposon insertion events

The pipeline is designed to be applied to model organisms, where reference annotation sets for SNPs/indels may not be available.
