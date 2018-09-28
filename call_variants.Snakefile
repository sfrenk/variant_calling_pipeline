import glob
import re
import sys
import os

############################    DESCRIPTION    ##############################

# Call Variants using GATK

##########################	NOTES	########################################

# Use setup_dir.sh -d <directory_name> -p gatk to set up the working directory

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = "/proj/ahmedlab/steve/seq/seq_data/trt1_polq1_rtel1_021218/"
EXTENSION = ".fastq.gz"
PAIRED = True
PROJECT_NAME = "cohort"

# The final list of variants must be filtered to remove false positives.
# There are two ways of doing this with C. elegans: 
# "hard" : use hard filters
# "bootstrap" : use highest scoring variants as a truth set for GATK's
#				machine learning-based filtering method.
# hard method works better when fewer samples (<10) are available
FITLER_METHOD = "hard"

# Remove items from the list below if you do not want the corresponding analysis to be performed.

# Options:
#	QC
#	snps_indels
#	structural_variants_svaba
#	structural_variants_meerkat (Useful for variants with single breakpoints eg. telomere recombination events)
#	CNV
#	repeat_copy_number
#	transposons
#	coverage
#	telomere_recombination
#	subtelomere_recombination

OPTIONS_LIST = ["QC", "coverage", "snps_indels", "structural_variants_svaba", "CNV", "repeat_copy_number", "transposons", "telomere_recombination", "subtelomere_recombination", "telomere_variants"]

############################	VARIABLES	################################

# Directory containing utilities for this pipeline
UTILS_DIR = "/nas/longleaf/home/sfrenk/pipelines/variant_calling"

#### Mapping and SNPs/indels ####

# Sequences of adapters to be trimmed
ADAPTERS="~/proj/seq/bbmap/adapters.fa"
# Fasta reference
REF = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.fa"
# BWA index
BWA_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bwa/genome.fa"
# Picard jar location
PICARD = "/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar"
# Picard reference dict for genome
REF_DICT = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.dict"
# snpEff directory
SNPEFF = "/nas/longleaf/home/sfrenk/local/src/snpEff/"

#### Structural Variants ####

# Breakpoint distance cutoff for defining unique events (note: this option is also used for CNVnator)
DISTANCE_CUTOFF=200

# Meerkat scripts directory
MEERKAT_SCRIPTS = "/nas/longleaf/home/sfrenk/local/Meerkat/scripts/"
# Meerkat genome directory
MEERKAT_GENOME = "/nas/longleaf/home/sfrenk/local/Meerkat/genomes/WS251"
# Directory containing samtools version 0.1.19
# Annoyingly, Meerkat requires this old version of samtools, but I want to keep the most recent samtools version in my PATH
SAMTOOLS_OLD = "/nas/longleaf/home/sfrenk/local/samtools/0.1.19/"

#### CNV ####
CNVNATOR = "~/local/CNVnator_v0.3/src/cnvnator"

#### Repeat Copy Number ####

# Tandem repeat loci gtf file
REPEATS_GTF = "~/proj/seq/WS251/repeats.gtf"

#### Transposons ####

# jitterbug directory
JITTERBUG_DIR = "/proj/ahmedlab/steve/Software/jitterbug/"
# GTF/GFF3 file containing all known transposon insertions in the genome
TRANSPOSONS_GTF = "/proj/ahmedlab/steve/seq/transposons/ce11_rebpase/ce11_transposons.gff3"
# Paired end support required
JITTERBUG_PAIRED = 2
# Soft-clipped support required
JITTERBUG_SOFT = 1

#### Telomere counts/ recombination ####

# Telomere sequence
TELOMERE_SEQ = "TTAGGC"
# Minimum number of telomere repeats required
TELOMERE_COUNT = 7
# Subtelomere bed file
SUBTELO_BED = UTILS_DIR + "/telomeres_2kb_subtelo.bed"

# motif_counter.sh location
MOTIF_COUNTER = "~/local/motif_counter.sh"

###############################################################################

SAMPLE_FILES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLE_FILES ]

# Paired end files must end in _1/_2 where 1 and 2 denote forward and reverse reads respectively. 
if PAIRED:
	SAMPLES = list(set([re.search("(^.+)_[12]$", x).group(1) for x in SAMPLES]))

# Check samples
if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

# Identify target output file(s)
output_files = {"QC" : "multiqc_report.html", "coverage" : expand("coverage/{sample}_coverage.txt.gz", sample = SAMPLES), "snps_indels" : "final/" + PROJECT_NAME + ".variants.ano.vcf", "structural_variants_svaba" : "svaba/" + PROJECT_NAME + ".log", "structural_variants_meerkat" : "meerkat/" + PROJECT_NAME + ".structural.meerkat.txt", "repeat_copy_number" : expand("repeats/{sample}_repeats.txt", sample = SAMPLES), "transposons" : "final/" + PROJECT_NAME + ".transposons.txt", "telomere_recombination" : "telomeres/telomere_mates.txt.temp", "subtelomere_recombination" : expand("subtelomeres/{sample}_summary.txt", sample = SAMPLES), "CNV" : "cnv/" + PROJECT_NAME + ".cnv.txt", "telomere_variants" : expand("telomere_variants/{sample}.txt", sample = SAMPLES)}


rule all:
	input:
		[ output_files[x] for x in OPTIONS_LIST ]

###############################################################################

#### Trimming and mapping reads ####

if PAIRED:

	rule trim:
		input:
			read1 = BASEDIR + "/{sample}_1" + EXTENSION,
			read2 = BASEDIR + "/{sample}_2" + EXTENSION
		output:
			out1 = "trimmed/{sample}_1.fastq",
			out2 = "trimmed/{sample}_2.fastq"
		params:
			adapter_file = ADAPTERS
		threads: 1
		log:
			"logs/{sample}_trim.log"
		shell:
			"module add bbmap; "
			"bbduk.sh -Xmx4g -ignorebadquality in1={input.read1} in2={input.read2} out1={output.out1} out2={output.out2} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

	rule fastqc:
		input:
			read1 = "trimmed/{sample}_1.fastq",
			read2 = "trimmed/{sample}_2.fastq"
		output:
			html1 = "metrics/fastq/{sample}_1_fastqc.html",
			html2 = "metrics/fastq/{sample}_2_fastqc.html",
			zip1 = "metrics/fastq/{sample}_1_fastqc.zip",
			zip2 = "metrics/fastq/{sample}_2_fastqc.zip"
		threads: 4
		log:
			"logs/{sample}_fastqc.log"
		shell:
			"module add fastqc; "
			"fastqc -o metrics/fastq {input.read1} {input.read2} &> {log}"

	rule bwa_mapping:
		input:
			trimmed1 = "trimmed/{sample}_1.fastq",
			trimmed2 = "trimmed/{sample}_2.fastq"
		output:
			"bwa_out/{sample}.sam"
		params:
			name = "{sample}",
			idx = BWA_INDEX
		log:
			"logs/{sample}_map.log"
		threads: 8
		log:
			"logs/{sample}_map.log"
		shell:
			"module add bwa; "
			"bwa mem -t {threads} -R '@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:ILLUMINA' {params.idx} {input.trimmed1} {input.trimmed2} > {output} 2> {log}"

else:
	rule trim:
		input:
			BASEDIR + "/{sample}" + EXTENSION
		output:
			"trimmed/{sample}.fastq"
		params:
			adapter_file = ADAPTERS
		threads: 1
		log:
			"logs/{sample}_trim.log"
		shell:
			"module add bbmap; "
			"bbduk.sh -Xmx4g -ignorebadquality in={input} out={output} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

	rule fastqc:
		input:
			"trimmed/{sample}.fastq"
		output:
			html = "metrics/fastq/{sample}_fastqc.html",
			zipfile = "metrics/fastq/{sample}_fastqc.zip"
		threads: 4
		log:
			"logs/{sample}_fastqc.log"
		shell:
			"module add fastqc; "
			"fastqc -o metrics/fastq {input} &> {log}"

	rule bwa_mapping:
		input:
			"trimmed/{sample}.fastq"
		output:
			"bwa_out/{sample}.sam"
		params:
			name = "{sample}",
			idx = BWA_INDEX
		log:
			"logs/{sample}_map.log"
		threads: 8
		shell:
			"module add bwa; "
			"bwa mem -t {threads} -R '@RG\\tID:{params.name}\\tSM:{params.name}\\tPL:ILLUMINA' {params.idx} {input} > {output} 2> {log}"

rule convert_to_bam:
	input:
		"bwa_out/{sample}.sam"
	output:
		"bam/{sample}.bam"
	params:
		name = "{sample}"
	log:
		"logs/{sample}_convert_to_bam.log"
	shell:
		#"sed -r 's/-R @RG\tID([^\t]+)\t/-R @RG\\t\1\\t/' {input} | samtools view -bh | samtools sort -o {output} -"
		"module add samtools; "
		"sed -r 's/-R @RG.*SM/-R @RG\tID:{params.name}\tSM/' {input} | sed 's/ID:bwa/ID:{params.name}/' | samtools view -bh - | samtools sort -o {output} - 2> {log}"
		#"samtools view -bh {input} | samtools sort -o {output} -"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai"
	shell:
		"module add samtools; "
		"samtools index {input}"

rule bam_metrics:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		"metrics/bam/{sample}/qualimapReport.html"
	params:
		output_dir="metrics/bam/{sample}"
	log:
		"logs/{sample}_bam_metrics.log"
	shell:
		"module add qualimap; "
		"qualimap bamqc -bam {input.bamfile} -outdir {params.output_dir} -gd HUMAN &> {log}"

#### Processing bam files ####

rule mark_duplicates:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		bam_out = "mrkdp/{sample}.bam",
		metrics = "metrics/picard/{sample}_metrics.txt"
	params:
		picard = PICARD
	log:
		"logs/{sample}_mark_duplicates.log"
	shell:
		"module add picard; "
		"java -jar {params.picard} MarkDuplicates INPUT={input.bamfile} OUTPUT={output.bam_out} METRICS_FILE={output.metrics} ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log}"

rule index_mrkdp_bam:
	input:
		"mrkdp/{sample}.bam"
	output:
		"mrkdp/{sample}.bam.bai"
	shell:
		"module add samtools; "
		"samtools index {input}"

#### QC compilation ###

if PAIRED:
	rule multiqc:
		input:
			fastq1 = expand("metrics/fastq/{sample}_1_fastqc.html", sample = SAMPLES),
			fastq2 = expand("metrics/fastq/{sample}_2_fastqc.html", sample = SAMPLES),
			bam = expand("metrics/bam/{sample}/qualimapReport.html", sample = SAMPLES),
			picard = expand("metrics/picard/{sample}_metrics.txt", sample = SAMPLES)
		output:
			"multiqc_report.html"
		log:
			"logs/multiqc.log"
		shell:
			"module purge; "
			"module add multiqc; "
			"multiqc -f metrics &> {log}"
else:
	rule multiqc:
		input:
			fastq = expand("metrics/fastq/{sample}_fastqc.html", sample = SAMPLES),
			bam = expand("metrics/bam/{sample}/qualimapReport.html", sample = SAMPLES),
			picard = expand("metrics/picard/{sample}_metrics.txt", sample = SAMPLES)
		output:
			"multiqc_report.html"
		log:
			"logs/multiqc.log"
		shell:
			"module purge; "
			"module add multiqc; "
			"multiqc -f metrics &> {log}"

#### Coverage ####

rule calculate_coverage:
	input:
		"mrkdp/{sample}.bam"
	output:
		"coverage/{sample}_coverage.txt.gz"
	log:
		"logs/{sample}_coverage.log"
	shell:
		"module add bedtools; "
		"bedtools genomecov -d -ibam {input} | gzip > {output} 2> {log}"


#### SNP/Indel calling ####

rule call_haplotypes_first_round:
	input:
		bamfile = "mrkdp/{sample}.bam",
		bamidx = "mrkdp/{sample}.bam.bai"
	output:
		"vcf_1/{sample}.g.vcf"
	params:
		ref = REF
	#threads: 8 Currently can't do multithreading for HaplotypeCaller
	log:
		"logs/{sample}_call_haplotypes1.log"
	shell:
		"module add gatk; "
		"gatk -T HaplotypeCaller -R {params.ref} -I {input.bamfile} --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -out_mode EMIT_ALL_CONFIDENT_SITES -o {output} 2> {log}"

# For samples > 200, need to combine gvcfs into batches with CombineGVCFs
# Have not yet properly built this into the pipeline

#rule combine_gvcfs_first_round:
#	input:
#		expand("vcf_1/{sample}.g.vcf", sample = SAMPLES)
#	output:
#		"vcf_1/" + PROJECT_NAME + ".first_pass.vcf"
#	params:
#		gvcf_list = GVCF_LIST1,
#		ref = REF
#	log:
#		"logs/combine_gvcfs1.log"
#	shell:
#		"gatk -T CombineGVCFs -R {params.ref} {params.gvcf_list} -o {output} 2> {log}"


## First round genotyping and base recalibration ##
# This will give an initial set of variants. The top variants will then be used for base quality score recalibration 

rule genotype_gvcfs_first_round:
	input:
		expand("vcf_1/{sample}.g.vcf", sample = SAMPLES)
	output:
		"vcf_1/" + PROJECT_NAME + ".first_pass.vcf"
	params:
		gvcf_list = " ".join([ "-V vcf_1/" + x + ".g.vcf" for x in SAMPLES ]),
		ref = REF
	log:
		"logs/genotype_gvcfs1.log"
	shell:
		"module add gatk; "
		"gatk -T GenotypeGVCFs -R {params.ref} {params.gvcf_list} -o {output} 2> {log}"

rule make_base_recalibration_table:
		input:
			bamfile = "mrkdp/{sample}.bam",
			gvcf = "vcf_1/" + PROJECT_NAME + ".first_pass.vcf"
			#recal_table = "recal/{sample}.recal.table"
		output:
			"recal/{sample}.recal.table"
		params:
			ref = REF
		log:
			"logs/{sample}_make_base_recalibration_table.log"
		shell:
			"module add gatk; "
			"gatk -T BaseRecalibrator -R {params.ref} -I {input.bamfile} --knownSites {input.gvcf} -o {output} 2> {log}"

rule base_recalibration:
	input:
		bamfile = "mrkdp/{sample}.bam",
		recal_table = "recal/{sample}.recal.table"
	output:
		"recal/{sample}.bqsr.bam"
	params:
		ref = REF
	log:
		"logs/{sample}_base_recalibration.log"
	shell:
		"module add gatk; "
		"gatk -T PrintReads -I {input.bamfile} -R {params.ref} -BQSR {input.recal_table} -o {output} 2> {log}"

rule index_recalibrated_bam:
	input:
		"recal/{sample}.bqsr.bam"
	output:
		"recal/{sample}.bqsr.bam.bai"
	shell:
		"module add samtools; "
		"samtools index {input}"

rule make_base_recalibration_after_table:
		input:
			bamfile = "recal/{sample}.bqsr.bam",
			gvcf = "vcf_1/" + PROJECT_NAME + ".first_pass.vcf"
			#recal_table = "recal/{sample}.recal.table"
		output:
			"recal/{sample}.after_recal.table"
		params:
			ref = REF
		log:
			"logs/{sample}_make_base_recalibration_after_table.log"
		shell:
			"module add gatk; "
			"gatk -T BaseRecalibrator -R {params.ref} -I {input.bamfile} --knownSites {input.gvcf} -o {output} 2> {log}"

rule make_recal_plots:
	input:
		before_table = "recal/{sample}.recal.table",
		after_table = "recal/{sample}.after_recal.table"
	output:
		"recal/{sample}.recal_plots.pdf"
	params:
		ref = REF
	log:
		"logs/{sample}_make_recal_plots.log"
	shell:
		"module add gatk; "
		"gatk -T AnalyzeCovariates -R {params.ref} -before {input.before_table} -after {input.after_table} -plots {output} 2> {log}"


## Second round genotyping ##

rule call_haplotypes_second_round:
	input:
		bamfile = "recal/{sample}.bqsr.bam",
		bamidx = "recal/{sample}.bqsr.bam.bai"
	output:
		"vcf_2/{sample}.g.vcf"
	params:
		ref = REF
	#threads: 8
	log:
		"logs/{sample}_call_haplotypes2.log"
	shell:
		"module add gatk; "
		"gatk -T HaplotypeCaller -R {params.ref} -I {input.bamfile} --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -out_mode EMIT_ALL_CONFIDENT_SITES -o {output} 2> {log}"

# For samples > 200, need to combine gvcfs into batches with CombineGVCFs
# Have not yet properly built this into the pipeline
#rule combine_gvcfs_second_round:
#	input:
#		expand("vcf_2/{sample}.g.vcf", sample = SAMPLES)
#	output:
#		"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
#	params:
#		gvcf_list = GVCF_LIST2,
#		ref = REF
#	log:
#		"logs/combine_gvcfs2.log"
#	shell:
#		"gatk -T CombineGVCFs -R {params.ref} {params.gvcf_list} -o {output} 2> {log}"

rule genotype_gvcfs_second_round:
	input:
		expand("vcf_2/{sample}.g.vcf", sample = SAMPLES)
	output:
		"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
	params:
		gvcf_list = " ".join([ "-V vcf_2/" + x + ".g.vcf" for x in SAMPLES ]),
		ref = REF
	log:
		"logs/genotype_gvcfs2.log"
	shell:
		"module add gatk; "
		"gatk -T GenotypeGVCFs -R {params.ref} {params.gvcf_list} -o {output} 2> {log}"


#### Variant filtering ####
# The variant set has been called, but must now be filtered to get rid of false positives. With human data, this is achieved using a validation set of known variants (eg dbSNP), which is used to train the machine learning algorithm in VariantRecalibrator. Such data doesn't exist for C elegans, so need to either set hard variables or bootstrap high-confidence variants form the set of called variants itself. In this pipeline, the top 10% highest quality variants are taken as the validation set for the bootstrap method.

if FITLER_METHOD == "hard":

	# Filter parameters are taken from the Cook et al. 2016 pot-2 Genetics paper

	rule extract_snps:
		input:
			"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
		output:
			"vcf_2/" + PROJECT_NAME + ".second_pass.snps.vcf"
		params:
			ref = REF
		log:
			"logs/extract_snps.log"
		shell:
			"module add gatk; "
			"gatk -T SelectVariants -R {params.ref} -V {input} -selectType SNP -o {output} &> {log}"


	rule filter_snps_hard:
		input:
			"vcf_2/" + PROJECT_NAME + ".second_pass.snps.vcf"
		output:
			"final/" + PROJECT_NAME + ".snps.vcf"
		params:
			ref = REF
		log:
			"logs/filter_snps_hard.log"
		shell:
			"module add gatk; "
			'gatk -T VariantFiltration -R {params.ref} -o {output} --variant {input} --filterExpression "DP <= 10 || MQ <= 40 || QUAL < 30 || DV/DP < 0.5" --filterName "hard_filter_snp" &> {log}'

	rule extract_indels:
		input:
			"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
		output:
			"vcf_2/" + PROJECT_NAME + ".second_pass.indels.vcf"
		params:
			ref = REF
		log:
			"logs/extract_indels.log"
		shell:
			"module add gatk; "
			"gatk -T SelectVariants -R {params.ref} -V {input} -selectType INDEL -o {output} &> {log}"


	rule filter_indels_hard:
		input:
			"vcf_2/" + PROJECT_NAME + ".second_pass.indels.vcf"
		output:
			"final/" + PROJECT_NAME + ".indels.vcf"
		params:
			ref = REF
		log:
			"logs/filter_indels_hard.log"
		shell:
			"module add gatk; "
			'gatk -T VariantFiltration -R {params.ref} -o {output} --variant {input} --filterExpression "DP <= 10 || MQ <= 40 || QUAL < 30 || DV/DP < 0.5" --filterName hard_filter_indel &> {log}' 

elif FITLER_METHOD == "bootstrap":
	rule get_snp_training_set:
		input:
			"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
		output:
			snps_training_file = "recal/training_snps.vcf",
		log:
			"logs/get_snp_training_set.log"
		params:
			utils_dir = UTILS_DIR
		shell:
			"module add vcftools picard; "
			"bash {params.utils_dir}/sample_vcf.sh -v {input} -o {output} 2> {log}"

	rule get_indel_training_set:
		input:
			"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
		output:
			snps_training_file = "recal/training_indels.vcf",
		log:
			"logs/get_indel_training_set.log"
		shell:
			"module add vcftools picard; "
			"sample_vcf -t indel -v {input} -o {output} 2> {log}"

	rule recal_snps:
		input:
			gvcf = "vcf_2/" + PROJECT_NAME + ".second_pass.vcf",
			training_set = "recal/training_snps.vcf"
		output:
			recal_file = "recal/snps.recal",
			tranches_file = "recal/snps.tranche",
			plots_file = "recal/snps.plots"
		params:
			ref = REF
		threads:
			8
		log:
			"logs/snp_recal.log"
		shell:
			"module add gatk; "
			"gatk -T VariantRecalibrator -nt {threads} -R {params.ref} -U ALLOW_SEQ_DICT_INCOMPATIBILITY --maxGaussians 4 -input {input.gvcf} -an QD -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode SNP -resource:highscoreset,known=true,training=true,truth=true,prior=10.0 {input.training_set} -recalFile {output.recal_file} -tranchesFile {output.tranches_file} -rscriptFile {output.plots_file} 2> {log}"

	rule recal_indels:
		input:
			gvcf = "vcf_2/" + PROJECT_NAME + ".second_pass.vcf",
			training_set = "recal/training_indels.vcf"
		output:
			recal_file = "recal/indels.recal",
			tranches_file = "recal/indels.tranche",
			plots_file = "recal/indels.plots"
		params:
			ref = REF
		threads:
			8
		log:
			"logs/indel_recal.log"
		shell:
			"module add gatk; "
			"gatk -T VariantRecalibrator -nt {threads} -R {params.ref} --maxGaussians 4 -U ALLOW_SEQ_DICT_INCOMPATIBILITY -input {input.gvcf} -an QD -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -resource:highscoreset,known=true,training=true,truth=true,prior=10.0 {input.training_set} -recalFile {output.recal_file} -tranchesFile {output.tranches_file} -rscriptFile {output.plots_file} 2> {log}"

	rule apply_snp_recal:
		input:
			gvcf = "vcf_2/" + PROJECT_NAME + ".second_pass.vcf",
			recal_file = "recal/snps.recal",
			tranches_file = "recal/snps.tranche"
		output:
			"final/" + PROJECT_NAME + ".snps.vcf"
		params:
			ref = REF
		log:
			"logs/apply_snp_recal.log"
		shell:
			"module add gatk; "
			"gatk -T ApplyRecalibration -R {params.ref} -input {input.gvcf} -tranchesFile {input.tranches_file} -recalFile {input.recal_file} -mode SNP --ts_filter_level 99.95 -o {output} 2> {log}"

	rule apply_indel_recal:
		input:
			gvcf = "vcf_2/" + PROJECT_NAME + ".second_pass.vcf",
			recal_file = "recal/indels.recal",
			tranches_file = "recal/indels.tranche"
		output:
			"final/" + PROJECT_NAME + ".indels.vcf"
		params:
			ref = REF
		log:
			"logs/apply_indel_recal.log"
		shell:
			"module add gatk; "
			"gatk -T ApplyRecalibration -R {params.ref} -input {input.gvcf} -tranchesFile {input.tranches_file} -recalFile {input.recal_file} -mode INDEL --ts_filter_level 95 -o {output} 2> {log}"

else:
	raise Exception("ERROR: FITLER_METHOD must be set to 'hard' or 'bootstrap'")

#### Annotating variants according to their predicted effect ####
# At this point, I also remove variants that didn't pass the filters. Amongst other things, this will allow the SNPs and indels to be merge at the end.

rule annotate_snps:
	input:
		"final/" + PROJECT_NAME + ".snps.vcf"
	output:
		"final/" + PROJECT_NAME + ".snps.ano.vcf"
	params:
		snpeff = SNPEFF
	threads:
		8
	log:
		"logs/annotate_snps.log"
	shell:
		"java -Xmx4g -jar {params.snpeff}/snpEff.jar -nodownload -t -c {params.snpeff}/snpEff.config WBcel235.86 {input} | grep -E '#|PASS' > {output} 2> {log}"

rule anotate_indels:
	input:
		"final/" + PROJECT_NAME + ".indels.vcf"
	output:
		"final/" + PROJECT_NAME + ".indels.ano.vcf"
	params:
		snpeff = SNPEFF
	threads:
		8
	log:
		"logs/annotate_snps.log"
	shell:
		"module add gatk; "
		"java -Xmx4g -jar {params.snpeff}/snpEff.jar -nodownload -t -c {params.snpeff}/snpEff.config WBcel235.86 {input} | grep -E '#|PASS' > {output} 2> {log}"

rule combine_variants:
	input:
		snps_vcf = "final/" + PROJECT_NAME + ".snps.ano.vcf",
		indels_vcf = "final/" + PROJECT_NAME + ".indels.ano.vcf"
	output:
		"final/" + PROJECT_NAME + ".variants.ano.vcf"
	params:
		ref = REF
	log:
		"logs/combine_variants.log"
	shell:
		"module add gatk; "
		"gatk -T CombineVariants -R {params.ref} --genotypemergeoption UNSORTED --variant {input.snps_vcf} --variant {input.indels_vcf} -o {output}"


#### Structural Variant calling ####

# I am currently using two different SV tools: meerkat and svaba

if "structural_variants_meerkat" in OPTIONS_LIST:
	rule symlink_bam:
		input:
			bamfile = "mrkdp/{sample}.bam",
			bamidx = "mrkdp/{sample}.bam.bai"
		output:
			bamfile = "meerkat/{sample}/{sample}.bam",
			bamidx = "meerkat/{sample}/{sample}.bam.bai"
		run:
			shell("ln -s ../../{input.bamfile} {output.bamfile}")
			shell("ln -s ../../{input.bamidx} {output.bamidx}")

	rule meerkat_pre_process:
		input:
			"meerkat/{sample}/{sample}.bam"
		output:
			preproc_log = "meerkat/{sample}/{sample}.pre.log",
			isinfo = "meerkat/{sample}/{sample}.isinfo"
		params:
			meerkat_scripts = MEERKAT_SCRIPTS,
			samtools_old = SAMTOOLS_OLD
		log:
			"logs/{sample}_meerkat_pre_process.log"
		shell:
			"module add perl; "
			"perl {params.meerkat_scripts}/pre_process.pl -S {params.samtools_old} -s 20 -k 1500 -q 15 -l 0 -b {input} 2> {log}"

	rule run_meerkat:
		input:
			bamfile = "meerkat/{sample}/{sample}.bam",
			preproc_log = "meerkat/{sample}/{sample}.pre.log",
			isinfo = "meerkat/{sample}/{sample}.isinfo"
		output:
			meerkat_bam = "meerkat/{sample}/{sample}.sr.bam",
			cluster_file = "meerkat/{sample}/{sample}.clusters"
		threads:
			8
		params:
			meerkat_scripts = MEERKAT_SCRIPTS,
			meerkat_genome = MEERKAT_GENOME,
			samtools_old = SAMTOOLS_OLD
		log:
			"logs/{sample}_run_meerkat.log"
		shell:
			"module add perl bwa; "
			"perl {params.meerkat_scripts}/meerkat.pl -S {params.samtools_old} -s 20 -d 5 -p 3 -o 1 -m 0 -l 0 -t {threads} -F {params.meerkat_genome}/fasta -b {input.bamfile} 2> {log}"

	rule meerkat_mechanism:
		input:
			bamfile = "meerkat/{sample}/{sample}.bam",
			meerkat_bam = "meerkat/{sample}/{sample}.sr.bam"
		output:
			"meerkat/{sample}/{sample}.variants"
		params:
			meerkat_scripts = MEERKAT_SCRIPTS,
			meerkat_genome = MEERKAT_GENOME
		log:
			"logs/{sample}_meerkat_mechanism.log"
		shell:
			"module add perl; "
			"perl {params.meerkat_scripts}/mechanism.pl -R {params.meerkat_genome}/*_rmsk.txt -b {input.bamfile} 2> {log}"

	# Combine meekat output into one VCF-style file
	rule compile_meerkat:
		input:
			expand("meerkat/{sample}/{sample}.variants", sample = SAMPLES)
		output:
			"meerkat/" + PROJECT_NAME + ".structural.meerkat.txt" 
		params:
			utils_dir = UTILS_DIR,
			distance_cutoff = DISTANCE_CUTOFF
		log:
			"logs/meerkat_compile.log"
		shell:
			"module add r; "
			"Rscript {params.utils_dir}/process_meerkat_output.R -d {params.distance_cutoff} -o {output} {input} 2> {log}"

if "structural_variants_svaba" in OPTIONS_LIST:
	rule run_svaba:
		input:
			bamfiles = expand("mrkdp/{sample}.bam", sample = SAMPLES),
			bamidxs = expand("mrkdp/{sample}.bam.bai", sample = SAMPLES)
		output:
			"svaba/" + PROJECT_NAME + ".log"
		params:
			ref = BWA_INDEX,
			input_list = " ".join([ "-t mrkdp/" + x + ".bam" for x in SAMPLES ]),
			output_base = "svaba/" + PROJECT_NAME
		threads:
			8
		log:
			"logs/svaba.log"
		shell:
			"module add bwa samtools; "
			"svaba run {params.input_list} -p {threads} -a {params.output_base} -G {params.ref} 2> {log}"

#### CNV ####
if "CNV" in OPTIONS_LIST:
	rule cnv_pre_process:
		input:
			bamfile = "mrkdp/{sample}.bam",
			bamidx = "mrkdp/{sample}.bam.bai"
		output:
			"cnv/{sample}"
		params:
			cnvnator = CNVNATOR,
			ref = REF
		log:
			"logs/{sample}_cnv_pre_process.log"
		shell:
			"{params.cnvnator} -root {output} -genome {params.ref} -tree {input.bamfile} &> {log} &&\
			{params.cnvnator} genome -root {output} -his 100 -d /nas/longleaf/home/sfrenk/proj/seq/WS251/genome/cnvnator/ &>> {log} &&\
			{params.cnvnator} -root {output} -stat 100 &>> {log} && \
			{params.cnvnator} -root {output} -partition 100 &>> {log}"

	rule cnv_call:
		input:
			"cnv/{sample}"
		output:
			"cnv/{sample}_cnvs.txt"
		params:
			cnvnator = CNVNATOR
		log:
			"logs/{sample}_cnv_call.log"
		shell:
			"{params.cnvnator} -root {input} -call 100 > {output} 2> {log}"

	rule cnv_compile:
		input:
			expand("cnv/{sample}_cnvs.txt", sample = SAMPLES)
		output:
			"cnv/" + PROJECT_NAME + ".cnv.txt" 
		params:
			utils_dir = UTILS_DIR,
			distance_cutoff = DISTANCE_CUTOFF
		log:
			"logs/cnv_compile.log"
		shell:
			"module add r; "
			"Rscript {params.utils_dir}/process_cnvnator_output.R -d {params.distance_cutoff} -o {output} {input} &> {log}"




#### Telomere/rDNA length analysis ####

# Note that pre-duplicate marked files are used here.

if "repeat_copy_number" in OPTIONS_LIST:
	rule repeat_copy_number:
		input:
			bamfile = "bam/{sample}.bam",
			bamidx = "bam/{sample}.bam.bai"
		output:
			"repeats/{sample}_repeats.txt"
		params:
			utils_dir = UTILS_DIR,
			repeats_gtf = REPEATS_GTF
		log:
			"logs/{sample}_copy_number.log"
		shell:
			"module add python samtools; "
			"python3 {params.utils_dir}/calculate_repeat_copy_number.py -g {params.repeats_gtf} -l -o {output} {input.bamfile} 2> {log}"

#### Transposons ####

if "transposons" in OPTIONS_LIST:
	rule find_transposon_insertions:
		input:
			bamfile = "mrkdp/{sample}.bam",
			bamidx = "mrkdp/{sample}.bam.bai"
		output:
			"transposons/{sample}/{sample}.TE_insertions_paired_clusters.gff3"
		params:
			jitterbug_dir = JITTERBUG_DIR,
			transposons_gtf = TRANSPOSONS_GTF,
			sample_name = "{sample}"
		threads:
			4
		log:
			"logs/{sample}_transposon_insertions.log"
		shell:
			"module load python/2.7.14 samtools; "
			"python {params.jitterbug_dir}/jitterbug.py --pre_filter -l {params.sample_name} -n {threads} -o transposons/{params.sample_name}/{params.sample_name} {input.bamfile} {params.transposons_gtf} &> {log}"

	rule jitterbug_filter:
		input:
			expand("transposons/{sample}/{sample}.TE_insertions_paired_clusters.gff3", sample = SAMPLES)
		output:
			"transposons/" + PROJECT_NAME + ".all.transposons.txt"
		params:
			utils_dir = UTILS_DIR,
			paired_support = JITTERBUG_PAIRED,
			soft_support = JITTERBUG_SOFT
		log:
			"logs/jitterbug_filter.log"
		shell:
			"module add python; "
			"python3 {params.utils_dir}/jitterbug_filter.py -p {params.paired_support} -s {params.soft_support} -o {output} {input} &> {log}"

	rule jitterbug_genotype:
		input:
			"transposons/" + PROJECT_NAME + ".all.transposons.txt"
		output:
			"final/" + PROJECT_NAME + ".transposons.txt"
		params:
			utils_dir = UTILS_DIR
		log:
			"logs/jitterbug_genotype.log"
		shell:
			"module add r; "
			"Rscript {params.utils_dir}/jitterbug_genotype.R -o {output} {input} &> {log}"


#### Telomere recombination ####

#if "telomere_recombination" in OPTIONS_LIST:
#	rule telomere_recombination:
#		input:
#			bamfile = "mrkdp/{sample}.bam",
#			bamidx = "mrkdp/{sample}.bam.bai"
#		output:
#			"telomere_recombination/{sample}.bam.bg"
#		params:
#			utils_dir = UTILS_DIR,
#			output_base = "telomere_recombination/{sample}.bam",
#			telomere_seq = TELOMERE_SEQ
#		log:
#			"logs/{sample}_telomere_recombination.log"
#		shell:
#			"python3 {params.utils_dir}/telomere_mates.py -p -o {params.output_base} -t {params.telomere_seq} -q 20 {input.bamfile} &> {log}"

#### Telomere count ####

if "telomere_recombination" in OPTIONS_LIST:
	rule telomere_counts:
		input:
			bamfiles = expand("mrkdp/{sample}.bam", sample = SAMPLES),
			bamidxs = expand("mrkdp/{sample}.bam.bai", sample = SAMPLES)
		output:
			counts_file = "telomeres.txt",
			bamfiles = expand("telomeres_BAM_FILES/telomeres_{sample}.bam_q0.bam", sample = SAMPLES),
		params:
			telomere_seq = TELOMERE_SEQ,
			count = TELOMERE_COUNT,
			motif_counter = MOTIF_COUNTER
		log:
			"logs/telomere_counts.log"
		shell:
			"module load samtools; "
			'printf "{params.telomere_seq}\n{params.count}\n" | bash {params.motif_counter} -i ./mrkdp -o telomeres -q 0 -Q 0 -p -u -v 2> {log} || true'

	rule sort_and_index_telomere_bam:
		input:
			"telomeres_BAM_FILES/telomeres_{sample}.bam_q0.bam"
		output:
			bamfile = "telomeres_BAM_FILES/telomeres_{sample}.sorted.bam",
			bamidx = "telomeres_BAM_FILES/telomeres_{sample}.sorted.bam.bai"
		log:
			"logs/{sample}_sort_and_index_telomere_bam.log"
		shell:
			"module add samtools; "
			"samtools sort -o {output.bamfile} {input} && samtools index {output.bamfile}"

	rule telomere_mates_bg:
		input:
			counts_file = "telomeres.txt",
			bamfiles = expand("telomeres_BAM_FILES/telomeres_{sample}.sorted.bam", sample = SAMPLES),
			bamidxs = expand("telomeres_BAM_FILES/telomeres_{sample}.sorted.bam.bai", sample = SAMPLES)
		output:
			"telomeres/telomere_mates.txt.temp"
		params:
			utils_dir = UTILS_DIR
		log:
			"logs/telomere_mates_bg.log"
		shell:
			"module load samtools bedtools; "
			"bash {params.utils_dir}/telomere_mates_bg.sh &> {log}"

if "subtelomere_recombination" in OPTIONS_LIST:
	rule subtelomere_mates:
		input:
			bamfile = "mrkdp/{sample}.bam",
			bamidx = "mrkdp/{sample}.bam.bai"
		output:
			"subtelomeres/{sample}_summary.txt"
		params:
			utils_dir = UTILS_DIR,
			output_base = "subtelomeres/{sample}",
			subtelo_bed = SUBTELO_BED
		log:
			"logs/{sample}_subtelomere_mates.log"
		shell:
			"module load samtools bedtools python; "
			"python3 {params.utils_dir}/subtelomere_mates.py -o {params.output_base} -b {params.subtelo_bed} {input.bamfile} &> {log}"

if "telomere_variants" in OPTIONS_LIST:
	rule telomere_variants:
		input:
			bamfile = "telomeres_BAM_FILES/telomeres_{sample}.sorted.bam",
			bamidx = "telomeres_BAM_FILES/telomeres_{sample}.sorted.bam.bai"
		output:
			"telomere_variants/{sample}.txt"
		params:
			utils_dir = UTILS_DIR
		log:
			"logs/{sample}_telomere_variants.log"
		shell:
			"module load python; "
			"python3 {params.utils_dir}/telomere_variants.py -o {output} {input.bamfile} &> {log}"
