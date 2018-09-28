import glob
import re
import sys
import os

############################    DESCRIPTION    ##############################

# Call Variants in human germline samples using GATK

##########################	NOTES	########################################

# Use setup_dir.sh -d <directory_name> -p gatk to set up the working directory

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ".fastq.gz"
PAIRED = True
PROJECT_NAME = "cohort"

# Remove items from the list below if you do not want the corresponding analysis to be performed.

# Options:
#	snps_indels
#	structural_variants
#	CNV
#	transposons
#	telomere_counts
#	telomere_variants

OPTIONS_LIST = ["snps_indels", "structural_variants", "CNV", "transposons", "telomere_counts", "telomere_variants"]

############################	VARIABLES	################################

# Directory containing utilities for this pipeline
UTILS_DIR = "/nas/longleaf/home/sfrenk/pipelines/variant_calling"

#### Mapping and SNPs/indels ####

# GATK module (gatk/{version})
GATK_MODULE="gatk/4.0.3.0"
# Sequences of adapters to be trimmed
ADAPTERS="~/proj/seq/bbmap/adapters.fa"
# Fasta reference
REF = "/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"
# BWA index
BWA_INDEX = "/proj/seq/data/hg38_UCSC/Sequence/BWAIndex/genome.fa"
# Picard jar location
PICARD = "/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar"
# Picard reference dict for genome
REF_DICT = "/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.dict"
# gatk bundle directory containing reference DBs for SNPs and indels
BUNDLE_DIR = "/nas/longleaf/home/sfrenk/proj/seq/human/hg38/gatk_bundle"
# snpEff directory
SNPEFF = "/nas/longleaf/home/sfrenk/local/src/snpEff/"

#### Structural Variants ####

#### Structural Variants ####
# svaba needs to be able to write a .fai file to the index directory
SVABA_BWA_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/human/hg38/bwa/genome.fa"

#### CNV ####
CNVNATOR = "~/local/CNVnator_v0.3/src/cnvnator"

# Breakpoint distance cutoff for defining unique events (note: this option is also used for CNVnator)
DISTANCE_CUTOFF=200

#### Repeat Copy Number ####

# Tandem repeat loci gtf file
REPEATS_GTF = "~/proj/seq/WS251/repeats.gtf"

#### Transposons ####

# jitterbug directory
JITTERBUG_DIR = "/proj/ahmedlab/steve/Software/jitterbug/"

# GTF/GFF3 file containing all known transposon insertions in the genome
TRANSPOSONS_GTF = "/proj/ahmedlab/steve/seq/human/hg38/transposons/hg38_transposons.gtf"

# Paired end support required
JITTERBUG_PAIRED = 2
# Soft-clipped support required
JITTERBUG_SOFT = 1

#### Telomere mates ####

# motif_counter.sh location
MOTIF_COUNTER = "~/local/motif_counter.sh"

# Telomere sequence
TELOMERE_SEQ = "TTAGGG"

# Minimum number of telomere repeats required
TELOMERE_COUNT = 6

# Subtelomere bed file
#SUBTELO_BED = UTILS_DIR + "/telomeres_2kb_subtelo.bed"

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
output_files = {"snps_indels" : "final/" + PROJECT_NAME + ".variants.ano.vcf", "structural_variants" : "svaba/" + PROJECT_NAME + ".log", "CNV" : "cnv/" + PROJECT_NAME + ".cnv.txt", "transposons" : "final/" + PROJECT_NAME + ".transposons.txt", "telomere_counts" : "telomeres.txt", "quality_check" : [expand("metrics/{sample}.insert_size_metrics", sample = SAMPLES), expand("metrics/{sample}.alignment_summary_metrics", sample = SAMPLES)], "telomere_variants" : expand("telomere_variants/{sample}.txt", sample = SAMPLES)}

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
			"module load bbmap; "
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
		log:
			"logs/{sample}_map.log"
		shell:
			"module load bwa; "
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
			"module load bbmap; "
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
		shell:
			"module load bwa; "
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
		"module load samtools; "
		"sed -r 's/-R @RG.*SM/-R @RG\tID:{params.name}\tSM/' {input} | sed 's/ID:bwa/ID:{params.name}/' | samtools view -bh - | samtools sort -o {output} - 2> {log}"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai"
	shell:
		"module load samtools; "
		"samtools index {input}"

#### bam quality metrics

rule get_bam_metrics:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		text_file = "metrics/{sample}.insert_size_metrics",
		pdf_file = "metrics/{sample}.insert_size_histogram.pdf"
	params:
		picard = PICARD,
		output_prefix = "metrics/{sample}"
	log:
		"logs/{sample}_get_bam_metrics.log"
	shell:
		"module load picard r; "
		'java -jar {params.picard} \
		      CollectMultipleMetrics \
		      INPUT={input.bamfile} \
		      OUTPUT={params.output_prefix} \
		      ASSUME_SORTED=true \
		      PROGRAM="null" \
		      PROGRAM="CollectBaseDistributionByCycle" \
		      PROGRAM="CollectInsertSizeMetrics" \
		      PROGRAM="MeanQualityByCycle" \
		      PROGRAM="QualityScoreDistribution" \
		      METRIC_ACCUMULATION_LEVEL="null" \
		      METRIC_ACCUMULATION_LEVEL="ALL_READS" 2> {log}'

rule get_alignment_summary:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		text_file = "metrics/{sample}.alignment_summary_metrics",
		pdf_file = "metrics/{sample}.gc_bias.pdf"
	params:
		picard = PICARD,
		ref = REF,
		output_prefix = "metrics/{sample}"
	log:
		"logs/{sample}_get_alignment_summary.log"
	shell:
		"module load picard r; "
		'java -jar {params.picard} \
			CollectMultipleMetrics \
			INPUT={input.bamfile} \
			REFERENCE_SEQUENCE={params.ref} \
			OUTPUT={params.output_prefix} \
			ASSUME_SORTED=true \
			PROGRAM="null" \
			PROGRAM="CollectAlignmentSummaryMetrics" \
			PROGRAM="CollectGcBiasMetrics" \
			METRIC_ACCUMULATION_LEVEL="null" \
			METRIC_ACCUMULATION_LEVEL="READ_GROUP" 2> {log}'

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
		"module load picard; "
		"java -jar {params.picard} MarkDuplicates INPUT={input.bamfile} OUTPUT={output.bam_out} METRICS_FILE={output.metrics} ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log}"

rule index_mrkdp_bam:
	input:
		"mrkdp/{sample}.bam"
	output:
		"mrkdp/{sample}.bam.bai"
	log:
		"logs/{sample}_index_mrkdp_bam.log"
	shell:
		"module load samtools; "
		"samtools index {input} &> {log}"


#### Coverage ####

rule calculate_coverage:
	input:
		"mrkdp/{sample}.bam"
	output:
		"coverage/{sample}_coverage.txt.gz"
	log:
		"logs/{sample}_coverage.log"
	shell:
		"module load bedtools; "
		"bedtools genomecov -d -ibam {input} | gzip > {output} 2> {log}"


#### Base quality score recalibration

rule make_base_recalibration_table:
		input:
			bamfile = "mrkdp/{sample}.bam"
		output:
			"recal/{sample}.recal.table"
		params:
			ref = REF,
			bundle_dir = BUNDLE_DIR,
			gatk_module = GATK_MODULE
		log:
			"logs/{sample}_make_base_recalibration_table.log"
		shell:
			"module add {params.gatk_module} python/2.7.12; "
			"gatk \
			BaseRecalibrator \
			-R {params.ref} \
			-I {input.bamfile} \
			--known-sites {params.bundle_dir}/hapmap_3.3.hg38.vcf.gz \
			--known-sites {params.bundle_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
			-O {output} &> {log}"

rule base_recalibration:
	input:
		bamfile = "mrkdp/{sample}.bam",
		recal_table = "recal/{sample}.recal.table"
	output:
		"recal/{sample}.bqsr.bam"
	params:
		gatk_module = GATK_MODULE,
		ref = REF
	log:
		"logs/{sample}_base_recalibration.log"
	shell:
		"module load {params.gatk_module} python/2.7.12; "
		"gatk \
		ApplyBQSR \
		-I {input.bamfile} \
		-R {params.ref} \
		-bqsr {input.recal_table} \
		-O {output} &> {log}"

rule index_recalibrated_bam:
	input:
		"recal/{sample}.bqsr.bam"
	output:
		"recal/{sample}.bqsr.bam.bai"
	log:
		"logs/{sample}_index_recalibrated_bam.log"
	shell:
		"module load samtools; "
		"samtools index {input}"


#### SNP/Indel calling ####

rule call_haplotypes:
	input:
		bamfile = "recal/{sample}.bqsr.bam",
		bamidx = "recal/{sample}.bqsr.bam.bai"
	output:
		"vcf/{sample}.g.vcf.gz"
	params:
		gatk_module = GATK_MODULE,
		ref = REF
	#threads: 8
	log:
		"logs/{sample}_call_haplotypes.log"
	shell:
		"module load {params.gatk_module} python/2.7.12; "
		"gatk \
		HaplotypeCaller \
		-R {params.ref} \
		-I {input.bamfile} \
		-O {output} \
		-ERC GVCF &> {log}"

# For samples > 200, need to combine gvcfs into batches with CombineGVCFs
# Have not yet properly built this into the pipeline
#rule combine_gvcfs_second_round:
#	input:
#		expand("vcf/{sample}.g.vcf", sample = SAMPLES)
#	output:
#		"vcf/" + PROJECT_NAME + ".second_pass.vcf"
#	params:
#		gvcf_list = GVCF_LIST2,
#		ref = REF
#	log:
#		"logs/combine_gvcfs2.log"
#	shell:
#		"gatk -T CombineGVCFs -R {params.ref} {params.gvcf_list} -o {output} 2> {log}"

# GATK4 can only genotype one file at a time, so need an extra step to merge gvcfs
rule merge_gvcfs:
	input:
		expand("vcf/{sample}.g.vcf.gz", sample = SAMPLES)
	output:
		"vcf/" + PROJECT_NAME + ".merged.gvcf.vcf.gz"
	params:
		gvcf_list = " ".join([ "--variant vcf/" + x + ".g.vcf.gz" for x in SAMPLES ]),
		gatk_module = GATK_MODULE,
		ref = REF
	log:
		"logs/merge_gvcfs.log"
	shell:
		"module load {params.gatk_module} python/2.7.12; "
		"gatk CombineGVCFs \
		--reference {params.ref} \
		{params.gvcf_list} \
		--output {output} \
		2> {log}"

rule genotype_gvcfs:
	input:
		"vcf/" + PROJECT_NAME + ".merged.gvcf.vcf.gz"
	output:
		"vcf/" + PROJECT_NAME + ".vcf.gz"
	params:
		gatk_module = GATK_MODULE,
		ref = REF
	log:
		"logs/genotype_gvcfs.log"
	shell:
		"module load {params.gatk_module} python/2.7.12 samtools; "
		"gatk \
		GenotypeGVCFs \
		--reference {params.ref} \
		--variant {input} \
		--output {output} 2> {log}"


#### Variant filtering ####
# The variant set has been called, but must now be filtered to get rid of false positives. With human data, this is achieved using a validation set of known variants (eg dbSNP), which is used to train the machine learning algorithm in VariantRecalibrator.

rule recal_snps:
	input:
		gvcf = "vcf/" + PROJECT_NAME + ".vcf.gz"
	output:
		recal_file = "recal/snps.recal",
		tranches_file = "recal/snps.tranche",
		plots_file = "recal/snps.plots"
	params:
		ref = REF,
		gatk_module = GATK_MODULE,
		bundle_dir = BUNDLE_DIR
	threads:
		2
	log:
		"logs/snp_recal.log"
	shell:
		"module load {params.gatk_module} r python/2.7.12; "
		"gatk \
		VariantRecalibrator \
		-R {params.ref} \
		-V {input.gvcf} \
		--resource hapmap,known=false,training=true,truth=true,prior=15.0:{params.bundle_dir}/hapmap_3.3.hg38.vcf.gz \
		--resource omni,known=false,training=true,truth=true,prior=12.0:{params.bundle_dir}/1000G_omni2.5.hg38.vcf.gz \
		--resource 1000G,known=false,training=true,truth=false,prior=10.0:{params.bundle_dir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
		--resource dbsnp,known=true,training=false,truth=false,prior=2.0:{params.bundle_dir}/dbsnp_146.hg38.vcf.gz \
		-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
		-mode SNP \
		-O {output.recal_file} \
		--tranches-file {output.tranches_file} \
		--rscript-file {output.plots_file} 2> {log}"

rule recal_indels:
	input:
		gvcf = "vcf/" + PROJECT_NAME + ".vcf.gz"
	output:
		recal_file = "recal/indels.recal",
		tranches_file = "recal/indels.tranche",
		plots_file = "recal/indels.plots"
	params:
		ref = REF,
		gatk_module = GATK_MODULE,
		bundle_dir = BUNDLE_DIR
	threads:
		2
	log:
		"logs/indel_recal.log"
	shell:
		"module load {params.gatk_module} r python/2.7.12; "
		"gatk \
		VariantRecalibrator \
		-R {params.ref} \
		--max-gaussians 4 \
		-V {input.gvcf} \
		--resource mills,known=false,training=true,truth=true,prior=12.0:{params.bundle_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
		--resource dbsnp,known=true,training=false,truth=false,prior=2.0:{params.bundle_dir}/dbsnp_146.hg38.vcf.gz \
		-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
		-mode INDEL \
		-O {output.recal_file} \
		--tranches-file {output.tranches_file} \
		--rscript-file {output.plots_file} 2> {log}"

rule apply_snp_recal:
	input:
		gvcf = "vcf/" + PROJECT_NAME + ".vcf.gz",
		recal_file = "recal/snps.recal",
		tranches_file = "recal/snps.tranche"
	output:
		"final/" + PROJECT_NAME + ".snps.vcf"
	params:
		gatk_module = GATK_MODULE,
		ref = REF
	log:
		"logs/apply_snp_recal.log"
	shell:
		"module load {params.gatk_module} python/2.7.12; "
		"gatk ApplyVQSR \
		-R {params.ref} \
		-V {input.gvcf} \
		--tranches-file {input.tranches_file} \
		--recal-file {input.recal_file} \
		-mode SNP \
		--truth-sensitivity-filter-level 99.5 \
		--create-output-variant-index true \
		-O {output} 2> {log}"

rule apply_indel_recal:
	input:
		gvcf = "vcf/" + PROJECT_NAME + ".vcf.gz",
		recal_file = "recal/indels.recal",
		tranches_file = "recal/indels.tranche"
	output:
		"final/" + PROJECT_NAME + ".indels.vcf"
	params:
		gatk_module = GATK_MODULE,
		ref = REF
	log:
		"logs/apply_indel_recal.log"
	shell:
		"module load {params.gatk_module} python/2.7.12; "
		"gatk ApplyVQSR \
		-R {params.ref} \
		-V {input.gvcf} \
		--tranches-file {input.tranches_file} \
		--recal-file {input.recal_file} \
		-mode INDEL \
		--truth-sensitivity-filter-level 99.0 \
		--create-output-variant-index true \
		-O {output} 2> {log}"

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
		"java -Xmx4g -jar {params.snpeff}/snpEff.jar -nodownload -t -c {params.snpeff}/snpEff.config hg38 {input} | grep -E '#|PASS' > {output} 2> {log}"

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
		"logs/annotate_indels.log"
	shell:
		"java -Xmx4g -jar {params.snpeff}/snpEff.jar -nodownload -t -c {params.snpeff}/snpEff.config hg38 {input} | grep -E '#|PASS' > {output} 2> {log}"

rule combine_variants:
	input:
		snps_vcf = "final/" + PROJECT_NAME + ".snps.ano.vcf",
		indels_vcf = "final/" + PROJECT_NAME + ".indels.ano.vcf"
	output:
		"final/" + PROJECT_NAME + ".variants.ano.vcf"
	params:
		picard = PICARD,
		ref = REF
	log:
		"logs/combine_variants.log"
	shell:
		"java -jar {params.picard} MergeVcfs \
          I={input.snps_vcf} \
          I={input.indels_vcf} \
          O={output}"


#### Structural Variant calling ####

if "structural_variants" in OPTIONS_LIST:
	rule run_svaba:
		input:
			bamfiles = expand("mrkdp/{sample}.bam", sample = SAMPLES),
			bamidxs = expand("mrkdp/{sample}.bam.bai", sample = SAMPLES)
		output:
			"svaba/" + PROJECT_NAME + ".log"
		params:
			ref = SVABA_BWA_INDEX,
			input_list = " ".join([ "-t mrkdp/" + x + ".bam" for x in SAMPLES ]),
			output_base = "svaba/" + PROJECT_NAME
		threads:
			8
		log:
			"logs/svaba.log"
		shell:
			"module load bwa samtools; "
			"svaba run {params.input_list} -p {threads} -a {params.output_base} -G {params.ref} 2> {log}"

#### CNV ####
if "CNV" in OPTIONS_LIST:

	rule make_cnvnator_genome_dir:
		input:
			REF
		output:
			"cnv_genome"
		params:
			utils_dir = UTILS_DIR
		log:
			"logs/make_cnvnator_genome_dir.log"
		shell:
			"module add python; "
			"python3 {params.utils_dir}/make_cnvnator_genome_dir.py -o {output} {input} &> {log}"


	rule cnv_pre_process:
		input:
			bamfile = "bam/{sample}.bam",
			bamidx = "bam/{sample}.bam.bai",
			cnv_genome_dir = "cnv_genome"
		output:
			"cnv/{sample}.root"
		params:
			cnvnator = CNVNATOR,
			ref = REF
		log:
			"logs/{sample}_cnv_pre_process.log"
		shell:
			"{params.cnvnator} -genome {params.ref} -root {output} -tree {input.bamfile} -unique; "
			"{params.cnvnator} -genome {params.ref} -root {output} -his 250 -d {input.cnv_genome_dir}; "
			"{params.cnvnator} -root {output} -stat 250; "
			"{params.cnvnator} -root {output} -partition 250 &> {log}"

	rule cnv_call:
		input:
			"cnv/{sample}.root"
		output:
			"cnv/{sample}_cnvs.txt"
		params:
			cnvnator = CNVNATOR
		log:
			"logs/{sample}_cnv_call.log"
		shell:
			"{params.cnvnator} -root {input} -call 250 > {output} 2> {log}"

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


#### Copy number analysis ####

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
			"module load python; "
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
			utils_dir = UTILS_DIR
		log:
			"logs/jitterbug_filter.log"
		shell:
			"module load python; "
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
			"module load r; "
			"Rscript {params.utils_dir}/jitterbug_genotype.R -o {output} {input} &> {log}"


#### Telomere count ####

if "telomere_counts" in OPTIONS_LIST:
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
