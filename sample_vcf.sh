#!/usr/bin/bash -e

# Variables
# picard jar location
picard="/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar"
# picard genome dictionary
ref_dict="/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.dict"

# Defaults
percentile=10
var_type=snp
outfile=sample.vcf
min_sample=0

usage="
	USAGE:
		Extract the top n% of variants of a VCF file by QUAL score
			-v/--vcf			Input VCF file
			-p/--percentile		Percentile to select (default: 10)
			-t/--type			Type of variant to choose (snp or indel. Default: snp)
			-o/--output			Output filename (default: sample.vcf)
			-m/--min_sample 	Minimum size of sample (eg. if 10% of variants is lower than this number, use this number instead. Default: 0)"


while [[ $# > 0 ]]
do
    key="$1"
    case $key in
    	-h|--help)
		echo "$usage"
		exit 0
		;;
    	-v|--vcf)
		vcf="$2"
		shift
		;;
        -p|--percentile)
		percentile="$2"
		shift
		;;
		-t|--type)
		var_type="$2"
		shift
		;;
		-o|--output)
		outfile="$2"
		shift
		;;
		-m|--min_sample)
		min_sample="$2"
		shift
		;;
    esac
shift
done

# Make temp directory for intermediate files

temp_dir=$(echo $outfile | sed -r 's/[/.]/_/g')

if [[ ! -d ${temp_dir} ]]; then
	mkdir ${temp_dir}
fi

# Select variant type

if [[ $var_type == "snp" ]]; then
	vcftools --remove-indels --vcf $vcf --recode --recode-INFO-all --stdout > ${temp_dir}/vars.vcf
elif [[ $var_type == "indel" ]]; then
	vcftools --keep-only-indels --vcf $vcf --recode --recode-INFO-all --stdout > ${temp_dir}/vars.vcf
else
	echo "ERROR: Invalid type. Please select snp or indel using -t/--type option"
fi

# Find the total number of variants
nvars=$(grep -v '#' ${temp_dir}/vars.vcf | awk '$6!="."' | wc -l)

# Calculate sample number
nsamp=$(awk "BEGIN { ns=(${percentile}/100)*${nvars}; i=int(ns); print (ns-i<0.5)?i:i+1 }")

if [[ nsamp -lt $min_sample ]]; then
	echo "Sample size ($nsamp) too small, sampling $min_sample variants instead."
	nsamp=$min_sample
fi

# Select the top n% variables
grep -v "#" ${temp_dir}/vars.vcf | awk '$6!="."' | sort -nrk 6,6 | head -$nsamp | sort -k 1,1 -k 2,2n > ${temp_dir}/vars_nsamp.vcf

# Add the original VCF header
grep "#" ${temp_dir}/vars.vcf | cat - ${temp_dir}/vars_nsamp.vcf > ${temp_dir}/unsorted.vcf

# Fix contig order
java -jar $picard SortVcf I=${temp_dir}/unsorted.vcf O=$outfile SEQUENCE_DICTIONARY=$ref_dict

# Remove temp directory
#rm -r ${temp_dir}
