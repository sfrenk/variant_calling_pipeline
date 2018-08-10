#!/usr/bin/bash

usage="
    Make bedgraph files for the discordant telomere read bam files created by motif_counter.sh

    ARGUMENTS
        -d/--dir
        Directory containing discordant telomere read bam files. (default: telomeres_BAM_FILES)

        -o/--output
        Output directory. This will be created if it doesn't already exist. (default: telomeres)

    "

# Defaults
dir="telomeres_BAM_FILES"
output_dir="telomeres"

# Arguments

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
    	-h|--help)
		echo "$usage"
		exit
		;;
        -d|--dir)
        dir="$2"
        shift
        ;;
        -o|--output)
		output_dir="$2"
		shift
		;;
    esac
	shift
done

if [[ ! -d $dir ]]; then
	echo "ERROR: $dir is not a directory. Please select a valid directory with the -d/--dir option"
	exit 1
fi

if [[ ! -d $output_dir ]]; then
	mkdir $output_dir
fi


shopt -s nullglob

files=(${dir}/*.bam)

for file in ${files[@]}; do
	
	base="$(echo $(basename $file) | sed -r 's/telomeres_(.+)\.bam_q0\.bam/\1/')"
	echo "Getting discordant telomere reads for ${base}"
	samtools sort -o ${output_dir}/${base}.bam $file
	samtools index ${output_dir}/${base}.bam
	bedtools genomecov -ibam ${output_dir}/${base}.bam -bg > ${output_dir}/${base}_telomere_mates.bg

done

touch ${output_dir}/telomere_mates.txt.temp	
