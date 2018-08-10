#!/usr/bin/env bash

# Hard variables

# Directory containing Snakemake and cluster.json files
snakedir='/nas/longleaf/home/sfrenk/pipelines/snakemake'

usage="Create directory with Snakemake files required for pipeline \n\n setup_dir -s <directory containing call_variants.Snakefile> -d <working directory (default: current directory)> \n\n"

dir="."

if [ -z "$1" ]; then
    echo "$usage"
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
        -h|--help)
		printf "$usage"
		exit
		;;
    esac
    shift
done


if [[ ! -d $dir ]]; then
	echo "ERROR: Invalid working directory"
	exit 1
fi

if [[ ! -d $snakefile_dir ]]; then
	echo "ERROR: Invalid snakefile directory"
	exit 1
fi

if [[ $dir == $snakefile_dir ]]; then
	echo "ERROR: Working directory must be in different directory to snakefile"
	exit 1
fi

# Copy over the snakefile
snakefile="call_variants.Snakefile"
cp ${snakefile_dir}/${snakefile} ./${snakefile}

# Edit base directory in Snakefile
# Remove trailing "/" from dir if it's there 
dir_name="$(echo $dir |sed -r 's/\/$//')"
sed -r -i -e "s,^BASEDIR.*,BASEDIR = \"${dir_name}\"," "$snakefile"

# Determine file extension
extension="$(ls $dir | grep -Eo "\.[^/]+(\.gz)?$" | sort | uniq)"

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

# Create Snakmake command script
printf "#!/usr/bin/bash\n" > "run_snakemake.sh"
printf "#SBATCH -t 2-0\n\n" >> "run_snakemake.sh"
printf "module add python\n\n" >> "run_snakemake.sh"
printf "snakemake -s $snakefile --rerun-incomplete --cluster-config ${snakedir}/cluster.json -j 100 --cluster \"sbatch -n {cluster.n} -N {cluster.N} -t {cluster.time}\"\n" >> "run_snakemake.sh"
