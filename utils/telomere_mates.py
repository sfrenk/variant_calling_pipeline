import argparse
import os
import pysam
import subprocess
import re
import sys

parser = argparse.ArgumentParser(description = "Create bam file containing non-telomeric mates of telomeric reads and provide telomere length estimate")

parser.add_argument("input_file", help = "input bam file")
parser.add_argument("-t", "--telomere_sequence", help = "Telomere repeat sequence (default:TTAGGC). Note that Ns are allowed", default = "TTAGGC")
parser.add_argument("-n", "--number_threshold", help = "Minimum number of telomere repeats required for a read to be classified as telomeric (default: 6)", type = int, default = 6)
parser.add_argument("-p", "--perfect", help = "Reads are only considered telomeric if they comprise a perfect run of repeats with no mismatches. -n option is ignored.", action = "store_true", default = False)
parser.add_argument("-o", "--output", help = "Output filename (default: telomere_discordant.bam)", default = "telomere_discordant.bam")
parser.add_argument("-q", "--min_mapq", help = "Minimum mapq for mate to be kept (default: 0)", type = int, default = 0)
parser.add_argument("-g", "--genome_size", help = "Size of organism reference genome in bases, used for normalizing telomere length estimate (default: 100286070 (C. elegans))", type = int, default = 100286070)


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


#### Telomere query construction ####

# "N" means any base, so convert to "." for regex
telo = args.telomere_sequence.replace("N", ".")

# Find reverse complement of telomere sequence
dna_dict = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
telo_revcomp = [ dna_dict[x] for x in reversed(telo) ]
telo_revcomp = "".join(telo_revcomp)
print("telomere_mates.py: searching for mates of reads consisting of perfect " + args.telomere_sequence + " or " + telo_revcomp.replace(".", "N") + " repeats")

def is_telomeric(read_seq = None, repeat_seq = None, perfect = False):
	# Function to determine whether read is telomeric

	if perfect == False:
		# Look for a given number of repeats and allow for gaps in between repeats
		result = re.match((".*" + telo) * args.number_threshold + ".*", read_seq)
	else:
		# Only accept perfect telomere sequence
		result = re.match(seq, telo_revcomp * 100)

	return(result)


#### Find telomeric reads ####

counter = 0

# Find telomeric reads in bam file
for read in bamfile.fetch():

	# See if read is telomeric
	seq = read.seq
	if is_telomeric(read.seq, telo, args.perfect) or is_telomeric(read.seq, telo_revcomp, args.perfect):
		
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


#### Telomere length estimation ####

print("telomere_mates.py: estimating total telomere length")

# Get total read count
flagstats = subprocess.check_output(["samtools", "flagstat", args.input_file])
total_read_count = int(re.search(b"([0-9]+) \+ [0-9]+ mapped", flagstats).group(1))

# Calculate estimate
# Note that genome size is divided by two, as we are counting fragments rather than reads

telomere_length_estimate = (counter * (args.genome_size/2)) / (len(args.telomere_sequence) * total_read_count)

# Write out results
est_file = re.sub("\.bam|\.bg", "", args.output)
with open(est_file, "w") as f:
	f.write(args.input_file + "\t" + str(telomere_length_estimate))

print("telomere_mates.py: done!")
