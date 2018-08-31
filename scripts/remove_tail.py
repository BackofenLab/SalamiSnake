#!/usr/bin/env python

import argparse
import logging
from sys import stdout
from Bio.SeqIO.QualityIO import FastqGeneralIterator
# avoid ugly python IOError when stdout output is piped into another program
# and then truncated (such as piping to head)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

tool_description = """
Remove a certain number of nucleotides from the 3'-tails of sequences in FASTQ
format.

Example usage:
- remove the last 7 nucleotides from file input.fastq, write result to file
  output.fastq:
remove_tail.py input.fastq 7 --out output.fastq
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "infile",
    help="Path to fastq file.")
parser.add_argument(
    "length",
    type=int,
    help="Remove this many nts.")
# optional arguments
parser.add_argument(
    "-o", "--outfile",
    help="Write results to this file.")
parser.add_argument(
    "-v", "--verbose",
    help="Be verbose.",
    action="store_true")
parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")

args = parser.parse_args()
if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
elif args.verbose:
    logging.basicConfig(level=logging.INFO, format="%(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
logging.info("  infile: '{}'".format(args.infile))
logging.info("  length: '{}'".format(args.length))
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

# check length parameter
if args.length < 0:
    raise ValueError("Length must be a positive integer, is '{}'.".format(args.length))

# remove tail
with (open(args.outfile, "w") if args.outfile is not None else stdout) as samout:
    for header, seq, qual in FastqGeneralIterator(open(args.infile)):

        # if removing tail would lead to an empty sequence,
        # set sequence to a single N to keep fastq synchronized
        if len(seq) <= args.length:
            seq = "N"
            qual = "B"
            logging.debug("read '{}' was too short to remove full tail".format(header))
            logging.debug("seq: {}".format(seq))
            logging.debug("len(seq): {}".format(len(seq)))
        else:
            seq = seq[0:-args.length]
            qual = qual[0:-args.length]

        samout.write("@%s\n%s\n+\n%s\n" % (header, seq, qual))
