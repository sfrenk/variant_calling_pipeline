#!/usr/bin/env bash

# Hard variables

# Directory containing Snakemake and cluster.json files
snakedir='/nas/longleaf/home/sfrenk/pipelines/snakemake'
dir=""

usage="\nCreate files required for running pipeline \n\n bash setup_dir.sh <options> \n\n\t-d/--dir directory containing fastq.gz files \n\t-s/--snakedir directory containing call_variants.Snakefile and utils directory \n\n You can then edit call_variants.Snakefile as required and run the pipeline: \n\n bash run_snakemake.\n\n"

if [ -z "$1" ]; then
    printf "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--dir)
        dir="$2"
        shift
        ;;
        -s|--snakedir)
		snakedir="$2"
		shift 
		;;
        -h|--help)
		printf "$usage"
		exit
		;;
    esac
    shift
done


if [[ ! -d $dir ]]; then
	echo "ERROR: Invalid directory"
	exit 1
fi

if [[ ! -d $snakedir ]]; then
	echo "ERROR: Invalid snakemake directory"
	exit 1
fi

if [[ $dir == "." ]] || [[ $dir == "./" ]]; then
	echo "ERROR: working directory must not be the same as snakemake directory"
	exit 1
fi

snakefile='call_variants.Snakefile'
modules="python bbmap bwa samtools gatk picard vcftools r perl bedtools"

# Copy over the snakefile
cp ${snakedir}/${snakefile} ./${snakefile}

# Edit base directory in Snakefile
base="$(basename ${dir})"
sed -r -i -e "s,^BASEDIR.*,BASEDIR = \"${dir}\"," "$snakefile"

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+" | sort | uniq)"

# Check if there are multiple file extensions in the same directory
ext_count="$(echo $extension | wc -l)"

if [[ ext_count == 0 ]]; then
	echo "ERROR: Directory is empty!"
elif [[ ext_count != 1 ]]; then
	echo "WARNING: Multiple file extensions found: using .fastq.gz"
	extension=".fastq.gz"
fi

# Edit extension in Snakefile
extension="\"${extension}\""
sed -i -r -e "s/^EXTENSION.*/EXTENSION = ${extension}/g" "$snakefile"
sed -i -r -e "s/^UTILS_DIR.*/UTILS_DIR = ${snakedir}\/utils/g" "$snakefile"

# Create Snakmake command script
printf "#!/usr/bin/bash\n\n" > "run_snakemake.sh"
printf "module add $modules\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --cluster-config ${snakedir}/utils/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> "run_snakemake.sh"
