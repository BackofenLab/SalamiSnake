import argparse
import pysam
import pandas
import logging

####################
##   ARGS INPUT   ##
####################

tool_description = """
The tool takes in a bed3 or interval file and returns
sequences in as fasta.
By default output is written to source file location.
Example usage:
fetch_DNA_sequence.py interval-file -o output.file
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "interval_file",
    help="Path to bed/interval file.")
parser.add_argument(
    "genome_file",
    help="Path to genome file (fasta).")
# optional arguments
parser.add_argument(
    "-o", "--outfile",
    help="Write results to this file.")
parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")

args = parser.parse_args()
if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
logging.info("  interval_file: '{}'".format(args.interval_file))
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

###################
##   READ DATA   ##
###################

print("[START]")
print("[NOTE] Read data")

# your interval data
file = args.interval_file
cl_regions = pandas.read_table(file, sep='\t', names=['chrom', 'start', 'stop'])

# link to the reference genome in fasta format
fastafile = pysam.Fastafile(args.genome_file)

print("[NOTE] finish")

#####################
##   WRITE OUTPUT  ##
#####################

print("[NOTE] Write output")

# check if outputfile was defined
outfile_name = ""
if args.outfile:
    outfile_name = args.outfile
else:
    outfile_name = "test-data/sequences.fa"
sequence_file = open(outfile_name, 'w')

# get sequence for coordinates
for i in range(0,len(cl_regions)):
    sequence_file.write(">" + str(cl_regions['chrom'][i]) + "_" + str(cl_regions['start'][i]) + "_" +  str(cl_regions['stop'][i]) + "\n")
    sequence_file.write(fastafile.fetch(str(cl_regions['chrom'][i]), int(cl_regions['start'][i]), int(cl_regions['stop'][i])) + '\n')

sequence_file.close()

print("[FINISH]")
