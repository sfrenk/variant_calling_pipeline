#!/usr/bin/bash

# Find snakefile and utils directories
script_dir="$(pwd)"
echo "Setting up variant calling pipelines in ${script_dir}"
snakefile_dir="${script_dir}/snakemake"
utils_dir="${snakefile_dir%/snakemake}/utils"

# Make setup_dir.sh executable
chmod +x ${utils_dir}/setup_dir.sh

# Edit setup_dir.sh to contain correct snakefile directory
sed -i -e "s|^snakefile_dir.*|snakefile_dir=${snakefile_dir}|g" ${utils_dir}/setup_dir.sh

# Edit Snakefiles to have the correct utils directory
snakefiles=("${snakefile_dir}/*.Snakefile")
utils_dir="\"${utils_dir}\""

for i in ${snakefiles[@]}; do
	sed -i -e "s|^UTILS_DIR.*|UTILS_DIR = ${utils_dir}|g" $i
done

echo "Done!"
