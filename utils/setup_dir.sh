#!/usr/bin/env bash

# Hard variables

# Directory containing Snakemake and cluster.json files
snakefile_dir='/nas/longleaf/home/sfrenk/pipelines/snakemake'

usage="Create directory with Snakemake files required for pipeline \n\n setup_dir -s <directory containing call_variants.Snakefile> -d <directory containing fastq files (default: current directory)> -p <pipeline (human or elegans)> \n\n"

dir="."
pipeline="elegans"

if [ -z "$1" ]; then
    printf "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -s|--snakefile_dir)
        snakefile_dir="$2"
        shift
        ;;
        -d|--dir)
        dir="$2"
        shift
        ;;
        -p|--pipeline)
        pipeline="$2"
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
	echo "ERROR: Invalid fastq directory"
	exit 1
fi

if [[ ! -d $snakefile_dir ]]; then
	echo "ERROR: Invalid snakefile directory"
	exit 1
fi

# Copy over the snakefile
if [[ $pipeline == "elegans" ]]; then
    snakefile="call_variants.Snakefile"
elif [[ $pipeline == "human" ]]; then
    snakefile="call_variants_human_germline.Snakefile"
else
    echo "ERROR: Invalid pipeline. Select human or elegans"
    exit 1
fi

cp ${snakefile_dir}/${snakefile} ./${snakefile}

# Edit base directory in Snakefile
# Remove trailing "/" from dir if it's there 
dir_name="$(echo $dir |sed -r 's/\/$//')"
dir_name=\"${dir_name}\"
sed -i -e "s|^BASEDIR.*|BASEDIR = ${dir_name}|g" "$snakefile"

# Determine file extension
# Check if there are multiple file extensions in the same directory
ext_count="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq | wc -l)"

if [[ $ext_count == 0 ]]; then
	echo "ERROR: Directory is empty!"
elif [[ $ext_count -gt 1 ]]; then
	echo "WARNING: Multiple file extensions found: using .fastq.gz"
	extension=".fastq.gz"
else
    extension="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq)"
fi

# Edit extension in Snakefile
extension="\"${extension}\""

sed -i -e "s|^EXTENSION.*|EXTENSION = ${extension}|g" "$snakefile"

# Create Snakmake command script
printf "#!/usr/bin/bash\n" > "run_snakemake.sh"
printf "#SBATCH -t 2-0\n\n" >> "run_snakemake.sh"
printf "module add python\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --rerun-incomplete --cluster-config ${snakefile_dir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> "run_snakemake.sh"
