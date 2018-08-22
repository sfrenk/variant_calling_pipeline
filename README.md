# Comprehensive variant calling pipeline for model organisms

This is a Snakemake pipeline for WGS data processing and analysis. The pipeline was designed to run on the LongLeaf cluster at UNC Chapel Hill, but can easily be configured to run on other clusters. Job submission and parallelization is handled by Snakemake. The following steps can be performed:

* fastq processing with bbduk.sh
* Read mapping with BWA-mem
* SNP/indel calling with GATK
* Functional annotation of SNPS/indels with snpEff
* Structural variant calling with SvABA and Meerkat
* Copy number estimation of tandem repeat elements
* CNV calling with CNVnator
* Identification of novel transposon insertion events with jitterbug
* Telomere copy number and recombination analysis with motif_counter.sh
* Telomere repeat variant analysis

The pipeline is designed to be applied to model organisms, where reference annotation sets for SNPs/indels may not be available.

## Dependencies

bbmap\
BWA\
samtools\
bedtools\
GATK v3.8\
Picard\
vcftools\
perl (Required for Meerkat structural variants)\
R\
Python v3\
Python v2 (Required for transposon insertion calling)\
Meerkat\
SvABA\
motif_counter.sh (https://sourceforge.net/projects/motifcounter/)

I currently use the following versions:

bbmap v38.12\
BWA v0.7.17\
samtools v1.8\
bedtools v2.26\
GATK v3.8-0\
Picard v2.10.3\
vcftools v0.1.15\
R v3.3.1\
perl v5.18.2\
python v3.5.1\
python v2.7.14\
Meerkat v0.189\
SvABA v0.2.1\
motif_counter.sh 2012-10-25 release\

## Usage

After downloading the script, check all the parameter variables in the call_variants.Snakefile (these are in capital letters and are found towards the beginning of the file). Whenever you run a pipeline, you can use the setup_dir.sh script in the utils directory to copy the snakefile into your working directory and create a SLURM-ready bash script for running the job (run_snakemake.sh).

```bash
# Set up working directory
bash setup_dir.sh -s <directory containing call_variants.Snakefile> -d <working directory (default: current directory)>

# Run pipeline
bash run_snakemake.sh

# Alternatively, send the run_snakemake.sh job itself to the cluster
sbatch run_snakemake.sh
```