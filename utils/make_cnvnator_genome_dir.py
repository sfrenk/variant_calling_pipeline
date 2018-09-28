#!/usr/bin/env python3
import re
import argparse
import sys
import os

parser = argparse.ArgumentParser(description = "Make genome directory (containing each chromosome sequence as a separate fasta file) for CNVnator")

parser.add_argument("input", help = "Reference genome fasta file")
parser.add_argument("-o", "--outdir", help = "output directory (default: current)", default = ".")

args = parser.parse_args()

# Make output directory if it doesn't already exist
if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

# Create placeholder for output file
outfile = open("null.txt", "w")

with open(args.input) as f:
	for line in f:
		if line.startswith(">"):
			#header
			header = line
			name = re.search(">(.+)", line.strip()).group(1)
			outfile.close()
			outfile = open(args.outdir + "/" + name + ".fa", "w")
			outfile.write(line)
		else:
			# Sequence
			outfile.write(line)

outfile.close()
