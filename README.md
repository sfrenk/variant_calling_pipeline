# Comprehensive variant calling pipeline for model organisms

This pipeline consists of:

* fastq processing with bbduk.sh
* Read mapping with BWA-mem
* SNP/indel calling with GATK
* Functional annotation of SNPS/indels with snpEff
* Structural variant calling with SvABA and Meerkat
* Copy number estimation of tandem repeat elements
* CNV calling with CNVnator
* Identification of novel transposon insertion events with jitterbug
* Telomere copy number and recombination analysis with motif_counter.sh

The pipeline is designed to be applied to model organisms, where reference annotation sets for SNPs/indels may not be available.

## Dependencies

bbmap
BWA
samtools
bedtools
GATK v3.8
Picard
vcftools
perl (Required for Meerkat structural variants)
R
Python v3
Python v2 (Required for transposon insertion calling)
motif_counter.sh (https://sourceforge.net/projects/motifcounter/)

I currently use the following versions:

bbmap v38.12
BWA v0.7.17
samtools v1.8
bedtools v2.26
GATK v3.8-0
Picard v2.10.3
vcftools v0.1.15
R v3.3.1
perl v5.18.2
python v3.5.1
python v2.7.14
motif_counter.sh 2012-10-25 release
