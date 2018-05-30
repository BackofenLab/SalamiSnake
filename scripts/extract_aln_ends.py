#!/usr/bin/env python

import argparse
import logging
from sys import stdout
from shutil import rmtree
from tempfile import mkdtemp
from pybedtools import BedTool
import pysam
# avoid ugly python IOError when stdout output is piped into another program
# and then truncated (such as piping to head)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

tool_description = """
Extract alignment ends from sam file and export to bed format.

The resulting bed file contains the outer coordinates of the alignments. The bed
name field is set to the read id and the score field is set to the edit distance
of the alignment. The crosslinked nucleotide is one nt upstream of the 5'-end of
the bed entries.

This script only reports results for alignments that are properly aligned in FR
("forward-reverse") direction.

By default output is written to stdout.

Input:
* alignments in SAM or BAM format (paired-end sequencing)

Output:
* bed6 file containing outer coordinates (sorted by read id)

Example usage:
- Extract coordinates from file input.bam and write to file output.bed
extract_aln_ends.py input.bam --out output.bed
"""


class DefaultsRawDescriptionHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                          argparse.RawDescriptionHelpFormatter):
    # To join the behaviour of RawDescriptionHelpFormatter with that of ArgumentDefaultsHelpFormatter
    pass


# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=DefaultsRawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "infile",
    help="Path to alignments in SAM or BAM format.")
# optional arguments
parser.add_argument(
    "-o", "--outfile",
    help="Write results to this file.")
# misc arguments
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
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

# convert to bam for use with pybedtools
try:
    # setup temporary directory
    tmpdir = mkdtemp()
    logging.debug("tmpdir: " + tmpdir)

    fn_sorted = tmpdir + "/sorted.bam"
    fn_fixedmates = tmpdir + "/fixedmates.bam"

    # sort by id
    logging.debug("calling samtools sort")
    pysam.sort(args.infile, "-n", "-o{}".format(fn_sorted), "-T{}".format(tmpdir))

    # fix mate information
    # also removes secondary and unmapped reads
    logging.debug("calling samtools fixmates")
    pysam.fixmate("-r", fn_sorted, fn_fixedmates)

    # bedtools bam2bed
    alns = BedTool(fn_fixedmates)
    alns_bedpe = alns.bam_to_bed(bedpe=True, mate1=True, ed=True)

    # determine alignment ends and write to file
    with (open(args.outfile, "w") if args.outfile is not None else stdout) as out:
        for i in alns_bedpe:
            chrom = i.fields[0]
            fmstart = i.fields[1]
            fmend = i.fields[2]
            smstart = i.fields[4]
            smend = i.fields[5]
            readid = i.fields[6]
            score = i.fields[7]
            fmstrand = i.fields[8]
            if fmstrand == "+":
                start = fmstart
                end = smend
            elif fmstrand == "-":
                start = smstart
                end = fmend
            else:
                logging.warning("Skipping {}, strand information is missing: '{}'".format(readid, i))
                continue
            out.write("\t".join([chrom, start, end, readid, score, fmstrand]) + "\n")
finally:
    logging.debug("removed tmpdir: " + tmpdir)
    rmtree(tmpdir)
