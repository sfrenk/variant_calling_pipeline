import glob
import re
import sys

############################    DESCRIPTION    ##############################

# Call Variants using GATK

# Required modules:
#	python bbmap bwa samtools gatk picard vcftools r perl bedtools


##########################	NOTES	########################################

# Use setup_dir.sh -d <directory_name> -p gatk to set up the working directory

# Base recalibration is currently not supported, so keep BASE_RECALIBRATION = False. I need to get a vcf file of known mutations in our lab strains.

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = "fastq2"
EXTENSION = ".fastq.gz"
PAIRED = True
PROJECT_NAME = "cohort"

# Remove items from the list below if you do not want the corresponding analysis to be performed
OPTIONS_LIST = ["coverage", "snps_indels", "structural_variants", "repeat_copy_number", "transposons", "telomere_recombination"]

############################	VARIABLES	################################

# Directory containing utilities for this pipeline
UTILS_DIR = "/nas/longleaf/home/sfrenk/pipelines/variant_calling"


#### Mapping and SNPs/indels ####

# Sequences of adapters to be trimmed
ADAPTERS="/nas/longleaf/apps/bbmap/37.62/bbmap/resources/adapters.fa"
# Fasta reference
REF = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.fa"
# BWA index
BWA_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bwa/genome"
# Picard jar location
PICARD = "/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar"
# Picard reference dict for genome
REF_DICT = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.dict"
# snpEff directory
SNPEFF = "/nas/longleaf/home/sfrenk/local/src/snpEff/"

#### Structural Variants ####

# Meerkat scripts directory
MEERKAT_SCRIPTS = "/nas/longleaf/home/sfrenk/local/Meerkat/scripts/"
# Meerkat genome directory
MEERKAT_GENOME = "/nas/longleaf/home/sfrenk/local/Meerkat/genomes/WS251"
# Directory containing samtools version 0.1.19
# Annoyingly, Meerkat requires this old version of samtools, but I want to keep the most recent samtools version in my PATH
SAMTOOLS_OLD = "/nas/longleaf/home/sfrenk/local/samtools/0.1.19/"

#### Copy Number ####

# Tandem repeat loci gtf file
REPEATS_GTF = "repeats.gtf"

#### Transposons ####

# jitterbug directory
JITTERBUG_DIR = "/proj/ahmedlab/steve/Software/jitterbug/"
# GTF/GFF3 file containing all known transposon insertions in the genome
TRANSPOSONS_GTF = "/proj/ahmedlab/steve/seq/transposons/ce11_rebpase/ce11_transposons.gff3"

#### Telomere mates ####

# Telomere coordinates bed file
TELOMERE_BED = "/nas/longleaf/home/sfrenk/proj/seq/telomere/telomeres.bed"

###############################################################################

SAMPLE_FILES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLE_FILES ]

# Paired end files must end in _1/_2 where 1 and 2 denote forward and reverse reads respectively. 
if PAIRED:
	SAMPLES = list(set([re.search("(^.+)_[12]$", x).group(1) for x in SAMPLES]))

# Construct input option for CombineGVCFs
GVCF_LIST1 = " ".join([ "-V vcf_1/" + x + ".g.vcf" for x in SAMPLES ])
GVCF_LIST2 = " ".join([ "-V vcf_2/" + x + ".g.vcf" for x in SAMPLES ])

if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

# Identify target output file(s)
output_files = {"coverage" : expand("coverage/{sample}_coverage.txt", sample = SAMPLES), "snps_indels" : "final/" + PROJECT_NAME + ".variants.ano.vcf", "structural_variants" : "final/" + PROJECT_NAME + ".structural.txt", "repeat_copy_number" : expand("repeats/{sample}_repeats.txt", sample = SAMPLES), "transposons" : "final/" + PROJECT_NAME + ".transposons.txt", "telomere_recombination" : expand("telomere_recombination/{sample}.bam.bg", sample = SAMPLES)}


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
			"bbduk.sh -Xmx4g -ignorebadquality in1={input.read1} in2={input.read2} out1={output.out1} out2={output.out2} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

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
		shell:
			"bwa mem -t {threads} -R '@RG\tID:{params.name}\tSM:{params.name}\tPL:ILLUMINA' {params.idx} {input.trimmed1} {input.trimmed2} > {output}"

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
			"bbduk.sh -Xmx4g -ignorebadquality in={input} out={output} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

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
		shell: "bwa mem -t {threads} -R '@RG\tID:{params.name}\tSM:{params.name}\tPL:ILLUMINA' {params.idx} {input} > {output} 2> {log}"

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
		"sed -r 's/-R @RG.*SM/-R @RG\tID:{params.name}\tSM/' {input} | sed 's/ID:bwa/ID:{params.name}/' | samtools view -bh - | samtools sort -o {output} - 2> {log}"
		#"samtools view -bh {input} | samtools sort -o {output} -"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai"
	shell:
		"samtools index {input}"


#### Processing bam files ####

rule mark_duplicates:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		bam_out = "mrkdp/{sample}.bam",
		metrics = "mrkdp/{sample}_metrics.txt"
	params:
		picard = PICARD
	log:
		"logs/{sample}_mark_duplicates.log"
	shell:
		"java -jar {params.picard} MarkDuplicates INPUT={input.bamfile} OUTPUT={output.bam_out} METRICS_FILE={output.metrics} ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log}"

rule index_mrkdp_bam:
	input:
		"mrkdp/{sample}.bam"
	output:
		"mrkdp/{sample}.bam.bai"
	shell:
		"samtools index {input}"


#### Coverage ####

rule calculate_coverage:
	input:
		"mrkdp/{sample}.bam"
	output:
		"coverage/{sample}_coverage.txt"
	log:
		"logs/{sample}_coverage.log"
	shell:
		"bedtools genomecov -d -ibam {input} > {output} 2> {log}"


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
		gvcf_list = GVCF_LIST1,
		ref = REF
	log:
		"logs/genotype_gvcfs1.log"
	shell:
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
		"gatk -T PrintReads -I {input.bamfile} -R {params.ref} -BQSR {input.recal_table} -o {output} 2> {log}"

rule index_recalibrated_bam:
	input:
		"recal/{sample}.bqsr.bam"
	output:
		"recal/{sample}.bqsr.bam.bai"
	shell:
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
		gvcf_list = GVCF_LIST2,
		ref = REF
	log:
		"logs/genotype_gvcfs2.log"
	shell:
		"gatk -T GenotypeGVCFs -R {params.ref} {params.gvcf_list} -o {output} 2> {log}"


#### Variant filtering ####
# The variant set has been called, but must now be filtered to get rid of false positives. With human data, this is achieved using a validation set of known variants (eg dbSNP), which is used to train the machine learning algorithm in VariantRecalibrator. Such data doesn't exist for C elegans, so need to bootstrap high-confidence variants form the set of called variants itself. In this pipeline, the top 10% highest quality variants are taken as the validation set.

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
		"bash {params.utils_dir}/sample_vcf.sh -v {input} -o {output} 2> {log}"

rule get_indel_training_set:
	input:
		"vcf_2/" + PROJECT_NAME + ".second_pass.vcf"
	output:
		snps_training_file = "recal/training_indels.vcf",
	log:
		"logs/get_indel_training_set.log"
	shell:
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
		"gatk -T ApplyRecalibration -R {params.ref} -input {input.gvcf} -tranchesFile {input.tranches_file} -recalFile {input.recal_file} -mode INDEL --ts_filter_level 95 -o {output} 2> {log}"

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
		"gatk -T CombineVariants -R {params.ref} --genotypemergeoption UNSORTED --variant {input.snps_vcf} --variant {input.indels_vcf} -o {output}"


#### Structural Variant calling ####

if "structural_variants" in OPTIONS_LIST:
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
			"perl {params.meerkat_scripts}/pre_process.pl -S {params.samtools_old} -s 20 -k 1500 -q 15 -l 0 -b {input} 2> {log}"

	rule run_meerkat:
		input:
			bamfile = "meerkat/{sample}/{sample}.bam",
			preproc_log = "meerkat/{sample}/{sample}.pre.log"
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
			"perl {params.meerkat_scripts}/mechanism.pl -R {params.meerkat_genome}/*_rmsk.txt -b {input.bamfile} 2> {log}"

	# Combine meekat output into one VCF-style file
	rule compile_meerkat:
		input:
			expand("meerkat/{sample}/{sample}.variants", sample = SAMPLES)
		output:
			"final/" + PROJECT_NAME + ".structural.txt" 
		params:
			meerkat_scripts = MEERKAT_SCRIPTS,
			samtools_old = SAMTOOLS_OLD,
			utils_dir = UTILS_DIR
		log:
			"logs/meerkat_compile.log"
		shell:
			"Rscript {params.utils_dir}/process_meerkat_output.R -o {output} {input} 2> {log}"

#### Telomere/rDNA length analysis ####

# Note that pre-duplicate marked files are used here.

if "tandem_repeat_copy_number" in OPTIONS_LIST:
	rule tandem_repeat_copy_number:
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
			"module load python/2.7.12; "
			"python {params.jitterbug_dir}/jitterbug.py --pre_filter -l {params.sample_name} -n {threads} -o transposons/{params.sample_name}/{params.sample_name} {input.bamfile} {params.transposons_gtf} 2> {log}"

	rule jitterbug_filter:
		input:
			expand("transposons/{sample}/{sample}.TE_insertions_paired_clusters.gff3", sample = SAMPLES)
		output:
			"transposons/" + PROJECT_NAME + ".all.transposons.txt"
		params:
			utils_dir = UTILS_DIR
		log:
			"logs/jitterbug_filter.log"
		shell:
			"python3 {params.utils_dir}/jitterbug_filter.py -o {output} {input} &> {log}"

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
			"Rscript {params.utils_dir}/jitterbug_genotype.R -o {output} {input} &> {log}"


#### Telomere recombination ####

if "telomere_recombination" in OPTIONS_LIST:
	rule telomere_recombination:
		input:
			bamfile = "mrkdp/{sample}.bam",
			bamidx = "mrkdp/{sample}.bam.bai"
		output:
			"telomere_recombination/{sample}.bam.bg"
		params:
			utils_dir = UTILS_DIR,
			output_base = "telomere_recombination/{sample}.bam",
			telomere_bed = TELOMERE_BED
		log:
			"logs/{sample}_telomere_recombination.log"
		shell:
			"python3 {params.utils_dir}/telomere_mates.py -o {params.output_base} -t {params.telomere_bed} -q 20 {input.bamfile} &> {log}"