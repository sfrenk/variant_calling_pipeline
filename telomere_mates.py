import argparse
import os
import pysam
import subprocess
import re
import sys

parser = argparse.ArgumentParser(description = "Create bam file containing non-telomeric mates of telomeric reads")

parser.add_argument("input_file", help = "input bam file")
parser.add_argument("-t", "--telomere_file", help = "Bed file containing telomere coordinates", required = True)
parser.add_argument("-o", "--output", help = "output filename (default: telomere_discordant.bam)", default = "telomere_discordant.bam")
parser.add_argument("-q", "--min_mapq", help = "minimum mapq for mate to be kept (default: 0)", type = int, default = 0)

args = parser.parse_args()

###############################################################################

## telomere coordinates ##

###############################################################################

#if os.getcwd() == bam_dir:
#	print("ERROR: Can't use bam directory as working directory!")
#	sys.exit()

###############################################################################

print("telomere_mates.py: starting...")

# Open input file
bamfile = pysam.AlignmentFile(args.input_file,"rb")

# Open output file
output_file = pysam.AlignmentFile(args.output,"wb", template = bamfile)

counter = 0

# Iterate through telomere coords to find telomere-mapped reads
with open(args.telomere_file) as f:

	# Get telomere coordinate
	for line in f.readlines():
		chrom = str(line.strip().split("\t")[0])
		start = int(line.strip().split("\t")[1])
		end = int(line.strip().split("\t")[2])

		for read in bamfile.fetch(chrom, start, end):

			# Find telomeric reads
			if read.is_paired:
				# Get mate, if possible
				try:
					mate = bamfile.mate(read)
				except ValueError:
					# Read is unmapped or has been filtered out
					continue

				# If the read pair is discordant, write the mate to the output file
				if not read.is_proper_pair and mate.mapq > args.min_mapq:
					output_file.write(mate)
					counter += 1

bamfile.close()
output_file.close()

print("telomere_mates.py: found " + str(counter) + " discordant telomere pairs. Writing to bg file...")
# Make bedGraph file of coverage for telomere mates

# First, need to sort and index bam file
subprocess.call(["samtools", "sort", "-o", args.output + ".sorted.bam", args.output])
subprocess.call(["samtools", "index", args.output + ".sorted.bam"])

# Make bedgraph
with open(args.output + ".bg", "w") as f:
	subprocess.call(["bedtools", "genomecov", "-ibam", args.output + ".sorted.bam", "-bg"], stdout = f)

print("telomere_mates.py: done!")
