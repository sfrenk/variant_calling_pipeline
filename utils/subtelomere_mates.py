import argparse
import os
import pysam
import subprocess
import re
import sys

parser = argparse.ArgumentParser(description = "Create bam file containing non-telomeric mates of telomeric reads")

parser.add_argument("input_file", help = "input bam file")
parser.add_argument("-b", "--bed", help = "Bed file containing subtelomere coordinates", required = True)
parser.add_argument("-o", "--output", help = "output file basename (default: subtelomere_discordant)", default = "subtelomere_discordant")
parser.add_argument("-q", "--min_mapq", help = "minimum mapq for subtelomeric read to be kept (default: 20)", type = int, default = 20)
parser.add_argument("-qm", "--min_mapq_mate", help = "minimum mapq for mate (default: 0)", type = int, default = 0)

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

# Keep track of the number of dischordant read pairs coming from different subtelomere regions

subtelo_counts = {}

# Iterate through subtelomere coords to find subtelomere-mapped reads
with open(args.bed) as f:

	# Get subtelomere coordinates
	for line in f.readlines():

		counter = 0

		chrom = str(line.strip().split("\t")[0])
		start = int(line.strip().split("\t")[1])
		end = int(line.strip().split("\t")[2])

		region = chrom + ":" + str(start) + "-" + str(end)

		print("Getting discordant reads for the region " + region)

		# Open output file
		output_file = pysam.AlignmentFile(args.output + "_" + region + ".bam","wb", template = bamfile)

		for read in bamfile.fetch(chrom, start, end):

			# Find subtelomeric reads
			if read.is_paired and read.mapq > args.min_mapq:
				# Get mate, if possible
				try:
					mate = bamfile.mate(read)
				except ValueError:
					# Read is unmapped or has been filtered out
					continue

				# If the read pair is discordant, write the mate to the output file
				if not read.is_proper_pair and mate.mapq > args.min_mapq_mate:
					output_file.write(mate)
					counter += 1

		output_file.close()
		subtelo_counts[region] = counter

bamfile.close()


print("telomere_mates.py: found the following number of discordant subtelomere pairs:\n\n")

with open(args.output + "_summary.txt", "w") as f:
	for x, count in subtelo_counts.items():
		print(x + "\t" + str(count))
		f.write(x + "\t" + str(count) + "\n")

# Process bam files and make bg files
print("\n\n")
for region in subtelo_counts.keys():
	print("processing:" + region)

	# Sort and index bam file
	subprocess.call(["samtools", "sort", "-o", args.output + "_" + region + ".sorted.bam", args.output + "_" + region + ".bam"])
	subprocess.call(["samtools", "index", args.output + "_" + region + ".sorted.bam"])

	# Make bedgraph
	with open(args.output + "_" + region + ".bg", "w") as f:
		subprocess.call(["bedtools", "genomecov", "-ibam", args.output + "_" + region + ".sorted.bam", "-bg"], stdout = f)

print("telomere_mates.py: done!")

