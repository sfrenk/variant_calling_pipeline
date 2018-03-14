#!/usr/bin/env python3

import argparse
import os
import re

parser = argparse.ArgumentParser(description = "Process the gff3 output files from jitterbug")

parser.add_argument("input", help = "Jitterbug output files to process", nargs = "+")
parser.add_argument("-o", "--output", help = "output file name name (default: jitterbug_output_cleaned.txt", default = "jitterbug_output_cleaned.txt")
parser.add_argument("-t", "--tag", help = "Identifier to be used in the fifth column (by default, this is extracted from gff3)", default = None)

args = parser.parse_args()

###############################################################################

output_file = open(args.output, "w")
output_gff = open(args.output+".gff", "w")

# Iterate through each input file

for file in args.input:

	gff = open(file)

	for line in gff:
		
		entry = line.strip().split("\t")

		# Get insertion coordinates
		chrom = entry[0]
		start = entry[3]
		end = entry[4]

		# Attirbutes are in the final column
		attributes = entry[8]

		# Get transposon names
		# Find origin of transposon anchor reads
		fwd = re.search("Inserted_TE_tags_fwd=([^;]*)", attributes).group(1).split(", ")
		rev = re.search("Inserted_TE_tags_rev=([^;]*)", attributes).group(1).split(", ")

		# Get number of supporting reads
		#fwd_support = re.search

		# If transposon insertion is legit, there should be both a forward and reverse tag with the same name. If more than one tags fits this criterium, the insertion is ambiguous and is not counted.
		
		paired_tags = []

		for tag in fwd:
			if tag in rev:
				paired_tags.append(tag)

		if len(paired_tags) == 1: #and len(fwd) == 1 and len(rev) == 1:

			# Define tag

			if not args.tag:
				tag = re.search("lib=([^;]+);", attributes).group(1)
			else:
				tag = args.tag

			output_file.write(chrom + "\t" + start + "\t" + end + "\t" + str(paired_tags[0]) + "\t" + tag + "\n")
			output_gff.write(line)

		else:

			# Group transposons into families. An insertion may be called assigned to multiple members of the same family due to sequence similarities

			fwd = [re.search("[A-Za-z]+", x).group() for x in fwd]
			rev = [re.search("[A-Za-z]+", x).group() for x in rev]

			for tag in fwd:
				if tag in rev:
					paired_tags.append(tag)

			if len(paired_tags) == 1: #and len(fwd) == 1 and len(rev) == 1:

			# Define tag

				if not args.tag:
					tag = re.search("lib=([^;]+);", attributes).group(1)
				else:
					tag = args.tag

				output_file.write(chrom + "\t" + start + "\t" + end + "\t" + str(paired_tags[0]) + "\t" + tag + "\n")
				output_gff.write(line)

	gff.close()
output_file.close()
output_gff.close()