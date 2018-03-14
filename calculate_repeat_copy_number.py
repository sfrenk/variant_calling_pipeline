import argparse
import os
import pysam
import subprocess
import re

parser = argparse.ArgumentParser(description = "Calculate coverage of regions from gtf file in a bam file")

parser.add_argument("bamfile", help = "bam file input")
parser.add_argument("-g", "--gff", help = "gtf file containing repeat regions")
parser.add_argument("-o", "--output", help = "Output filename (default: repeat_counts.txt)", default = "repeat_counts.txt")
parser.add_argument("-l", "--length", help = "By default, repeat length is calculated by subtracting the start from the end coordinate. Use this option if the length of the repeat is supplied as an attribute in the gff file", action = "store_true", default = False)
parser.add_argument("-s", "--genome_size", help = "size of organism reference genome in bases (default: 100286070(c elegans))", type = int, default = 100286070)

args = parser.parse_args()

###############################################################################
## VARIABLES ##

# Minimum mapping quality
# NOTE: if this is changed, then need to change the total read count.

min_mapq = 0

# Attribute string. This is used to select the attribute containing the feature name. This must contain any white space preceding the attribute name, but quotes should be excluded. Eg: for gene_id "pot-1", att string = "gene_id "

att_string = 'gene_id "?([^ \t;"]+)'


###############################################################################

output_file = open(args.output, "w")

# Get a list of region names
# The gff file contains coordinates for all loci to be analyzed

gff = open(args.gff, "r")

###############################################################################

print("Processing " + str(args.bamfile))

# Get total read count

flagstats = subprocess.check_output(["samtools", "flagstat", args.bamfile])

total_read_count = int(re.search(b"([0-9]+) \+ [0-9]+ mapped", flagstats).group(1))
#paired_pairs = int(re.search(b"([0-9]+) \+ [0-9]+ with itself and mate mapped", flagstats).group(1))
#total_read_count = total_reads - (paired_pairs/2)

# Open the bam file

bamfile_reads = pysam.AlignmentFile(args.bamfile,"rb")

# Create dictionaries for length and read count. These will be populated by the regions in the gff file

lengths = {}
counts = {}

for line in gff:

	# Define region
	
	chromosome = line.strip().split("\t")[0]
	start = int(line.strip().split("\t")[3])
	end = int(line.strip().split("\t")[4])
	
	try:
		name = re.search(att_string, line).group(1)
	except AttributeError:
		name = str(chromosome) + "_" + str(start) + "_" + str(end)

	print("Getting copy number for " + name)

	# Enter region into dictionaries if it isn't already there

	if name not in lengths:
		lengths[name] = 0
		counts[name] = 0

	
	# Count reads that map within the gff region

	read_names = []
	read_counts = 0

	read_names = [read.qname for read in bamfile_reads.fetch(chromosome, start, end) if read.mapq >= min_mapq]
	read_counts = len(read_names)

	counts[name] += read_counts

	# Record region length

	if args.length:
		region_length = int(re.search("length[=\s]([0-9]+)", line).group(1))
		lengths[name] = region_length
	else:
		region_length = end - start	
		lengths[name] += region_length

bamfile_reads.close()


# Calculate region copy number
	# Use the formula N = RG/LT
		# R = number of read pairs in region
		# G = total genome size in bp (100286070)
		# L = size of region
		# T = total number of read pairs in project

# cn = (read_counts * 100286070) / (region_length * total_read_count)

copy_numbers = {}

for k, v in counts.items():
	copy_numbers[k] = (counts[k] * args.genome_size) / (lengths[k] * total_read_count)

# Write out copy number results for the strain to the output file
for k, v in copy_numbers.items():
	output_file.write(str(k) + "\t" + str(v) + "\n")

gff.close()
output_file.close()
