
import argparse
import subprocess as sb
import matplotlib.pyplot as plt
import os
import numpy
import math
import fit_nbinom as fnb
import sys

#########################
##   NECESSARY TOOLS   ##
#########################

# bedtools

####################
##   ARGS INPUT   ##
####################

tool_description = """
The tool can evalute the profile of peaks.
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description, usage='%(prog)s [options]',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

# version
parser.add_argument(
    "-v", "--version", action="version", version="%(prog)s 0.0")

# positional arguments
parser.add_argument(
    "-a", "--input_bed",
    metavar='*.bed',
    required=True,
    help="Paths to the peak file in bed6 format.")
parser.add_argument(
    "-b", "--input_bam",
    metavar='*.bai',
    required=True,
    help="List of paths to the read bam files used for the peakcalling.")

# optional arguments
parser.add_argument(
    "-o", "--output_folder",
    metavar='path/',
    default=os.getcwd(),
    help="Write results to this path.")

######################
##   CHECKS INPUT   ##
######################

parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")

args = parser.parse_args()

# Check if peak file is in bed6 format
bedfile = open(args.input_bed, "r")
firstline = bedfile.readline()
if ( len(firstline.strip("\n").split("\t")) <  6 ):
    sys.exit("[ERROR] Peakfile has to be in bed6 format!")
bedfile.close()

##############
##   FUNC   ##
##############

def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

##############
##   MAIN   ##
##############

outfilename = args.input_bam.split("/")
outfilename = outfilename[len(outfilename)-1]
outfilename = outfilename.strip(".bam")
outfilename = outfilename.strip(".bed")

# Generate Coverage file
coverage_file_name = "{}/{}_coverage.tsv".format(args.output_folder, outfilename)
sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(args.input_bed, args.input_bam, coverage_file_name), shell=True).wait()

peaks_file = open(args.input_bed, "r")
num_peaks = get_line_count(peaks_file)
peaks_file.close()

print("[NOTE] {} peaks will be evaluated.".format(num_peaks))

mean_coverage_peaks_dict = dict()
prob_success_peaks_dict = dict()
variance_coverage_peaks_dict = dict()
num_bp_peaks_dict = dict()
coordinates_dict = dict()

# Evaluate peak profiles for each replicate

coverage_file = open(coverage_file_name, "r")
num_coverage_lines = get_line_count(coverage_file)
coverage_file.close()

coverage_file = open(coverage_file_name, "r")

peak_counter = -1
peak_cov_list = []

# Calcualte mean and variance of peak coverage profiles
line_count = 0
for line in coverage_file:
    line_count += 1
    data = line.strip("\n").split("\t")
    bp = int(data[len(data)-2])
    cov = int(data[len(data)-1])

    if(bp == 1):
        if( peak_counter != -1 ):
            #mean_coverage_peaks[peak_counter] = numpy.mean(peak_cov_list)
            #variance_coverage_peaks[peak_counter] = numpy.std(peak_cov_list)
            # The fit is the alternative version of the NB. But I get the expected number of successes and the
            # probability of success.
            if ( not all(v == 0 for v in peak_cov_list) ):
                nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                mean_coverage_peaks_dict[peak_counter] = nb_fit["size"]
                prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1-nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
            else:
                mean_coverage_peaks_dict[peak_counter] = 0.0
                prob_success_peaks_dict[peak_counter] = 0.0
                variance_coverage_peaks_dict[peak_counter] = 0.0
        peak_cov_list = []
        peak_cov_list.append(cov)
        peak_counter += 1
        coordinates_dict[peak_counter] = [data[0], data[1], data[2]]
    else:
        peak_cov_list.append(cov)
        num_bp_peaks_dict[peak_counter] = bp

        if ( line_count == num_coverage_lines ):
            if (  not all(v == 0 for v in peak_cov_list) ):
                nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                mean_coverage_peaks_dict[peak_counter] = nb_fit["size"]
                prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1 - nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
            else:
                mean_coverage_peaks_dict[peak_counter] = 0.0
                prob_success_peaks_dict[peak_counter] = 0.0
                variance_coverage_peaks_dict[peak_counter] = 0.0

coverage_file.close()

filtered_num_peaks = 0
varcoeff_coverage_peaks_dict = dict()

# Calcualte Variantioncoefficient of peak coverage profile
for i in range(0, num_peaks):

    if (mean_coverage_peaks_dict[i] > 0):
        varcoef = 1 / math.sqrt(mean_coverage_peaks_dict[i] * (1 - prob_success_peaks_dict[i]))

        if ( math.isnan(varcoef) ):
            print(varcoef)

        norm_varvoef = varcoef / math.sqrt(num_bp_peaks_dict[i]-1) # Taking the estimation of the standard deviation into account
        varcoeff_coverage_peaks_dict[i] = norm_varvoef

        filtered_num_peaks += 1
    else:
        varcoeff_coverage_peaks_dict[i] = -0.01

print("[NOTE] {} peaks are covered.".format(filtered_num_peaks))

# Filter our peaks that are completly uncovered
filtered_varcoeff_coverage_peaks = []
for i in varcoeff_coverage_peaks_dict.values():
    if (i >= 0):
        filtered_varcoeff_coverage_peaks.append(i)


# Normalize all VC so that the scale goes from 0 to 1 only
one = numpy.max(filtered_varcoeff_coverage_peaks)

for i in range(0, len(filtered_varcoeff_coverage_peaks)):
    filtered_varcoeff_coverage_peaks[i] = filtered_varcoeff_coverage_peaks[i]/one

print("[NOTE] Generate Plot")

# Make vase plot of variationkoefficients
f = plt.figure()
plt.violinplot(filtered_varcoeff_coverage_peaks)
plt.ylabel('Normalized Variationcoefficient of the Peak Profiles')
f.savefig(args.output_folder + "/VC_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')

print("[NOTE] Generate Output Tabular")
out_tab_file_name = args.output_folder + "/VC_tab_{}_tmp.bed".format(outfilename)
out_tab_file = open(out_tab_file_name, "w")
for i in range(0, num_peaks):
    coords = coordinates_dict[i]
    out_tab_file.write("{}\t{}\t{}\tpeak_{}\t{}\t.\t{}\t{}\n".format(coords[0], coords[1], coords[2],
                                                             i, varcoeff_coverage_peaks_dict[i],
                                                             num_bp_peaks_dict[i], mean_coverage_peaks_dict[i]))
out_tab_file.close()

sb.Popen("sort -r -k 5 {} > {}".format(out_tab_file_name, "{}/VC_tab_{}.bed".format(args.output_folder, outfilename)), shell=True).wait()
sb.Popen("rm {}".format(out_tab_file_name), shell=True).wait()

# testset = numpy.random.negative_binomial(10, .8, 10000)
# f = plt.figure()
# plt.hist(testset)
# f.savefig(args.output_folder + "/test1.pdf", bbox_inches='tight')
#
# print(fnb.fit_nbinom(testset))