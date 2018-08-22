#!/usr/bin/env python3

import os
import re
import pysam
import argparse
import sys

def make_telomere_variant_dictionary(telo_repeat):

	'''
	Make dictionary for every possible telomere variant combination
	'''

	telo_repeat = telo_repeat.upper()

	if re.match("[^ATCG]", telo_repeat):
		print("ERROR: telomere repeat should only contain the following characters: ATCG")
		sys.exit(1)

	base_dict = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}

	# Make table for telomere variant counts
	variant_counts_dict = {}
	for i in range(len(telo_repeat)):
		
		# Get each combination of telomere repeat in which one base is modified
		for j in base_dict.keys():
			var = [ x for x in telo_repeat ]
			var[i] = j

			# also get reverse complement
			var_revcomp = [ base_dict[x] for x in var ][::-1]

			# Add repeats to dictionary
			var = "".join(var)
			variant_counts_dict[var] = 0
			
			var_revcomp = "".join(var_revcomp)	
			variant_counts_dict[var_revcomp] = 0
	return(variant_counts_dict)

def iterate_bam(filename, variant_counts_dict):
	
	'''
	Count the number of occurrences of each repeat variant in a bam file
	'''

	# Open input file
	bamfile = pysam.AlignmentFile(filename,"rb")

	for read in bamfile.fetch():
		for telo_rep in variant_counts_dict.keys():
			matches = re.findall(telo_rep, read.seq.upper())
			variant_counts_dict[telo_rep] = variant_counts_dict[telo_rep] + len(matches)

	return(variant_counts_dict)

	bamfile.close()


def main():

	parser = argparse.ArgumentParser(description = "Get counts for the number of occurrences of telomere variants in a bam file. This script will look for every possible sequence combination in which one base is changed (including reverse complement.")

	parser.add_argument("input_file", help = "input bam file")
	parser.add_argument("-t", "--telomere_repeat", help = "Telomere repeat sequence (default: TTAGGC)", default = "TTAGGC")
	parser.add_argument("-o", "--output", help = "output filename (default: telomere_variant_counts.txt)", default = "telomere_variant_counts.txt")

	args = parser.parse_args()

	telomere_repeat_counts = make_telomere_variant_dictionary(args.telomere_repeat)

	telomere_repeat_counts = iterate_bam(args.input_file, telomere_repeat_counts)

	with open(args.output, "w") as f:
		for key, val in telomere_repeat_counts.items():
			f.write("\t".join([key, str(val)]) + "\n")

if __name__ == "__main__":
	main()
	