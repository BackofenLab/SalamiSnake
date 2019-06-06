
import argparse
import subprocess as sb
import matplotlib.pyplot as plt
import os
import numpy
import math
import fit_nbinom as fnb
import sys

plt.switch_backend('agg')

#########################
##   NECESSARY TOOLS   ##
#########################

# bedtools

##############
##   FUNC   ##
##############

# Function to obtain the number of lines in a file.
def get_line_count(file):
    count = 0
    for line in file:
        count += 1
    return count

def main():

    ####################
    ##   ARGS INPUT   ##
    ####################

    tool_description = """
    The tool can evalute the profile of peaks. Provide as an input the peaks you want to evalutate
    in bed6 format and the reads you used for the peak detection in bed or bam format. The user obtains
    a distributions of the variationcoefficient (VC) which can be used to evaluate the profile landscape. 
    In addition, the tool generates a list ranking the peaks based on the VC. The VC range from 0.0 for a broad 
    to 1.0 for a sharp peak. The ranked list has the following columns: chr, start, end, peakid, VC, strand,
    peak length, mean read coverage. 
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
        help="Path to the peak file in bed6 format.")
    parser.add_argument(
        "-b", "--input_bam",
        metavar='*.bai',
        required=True,
        help="Path to the read bam file used for the peakcalling in bed or bam format.")
    parser.add_argument(
        "-c", "--chr_file",
        metavar='*.txt',
        required=True,
        help="Path to the chromosome length file.")

    # optional arguments
    parser.add_argument(
        "-o", "--output_folder",
        metavar='path/',
        default=os.getcwd(),
        help="Write results to this path.")
    parser.add_argument(
        "--length_norm",
        action='store_true',
        help="Set length normalization. StoatyDive will expand every peak to the maximal length.")
    parser.add_argument(
        "--length_norm_value",
        metavar='int',
        help="Set length normalizatiuon value (maximum peak length).")
    parser.add_argument(
        "--max_norm_value",
        metavar='float',
        help="Provide a maximum value for the coefficient of variation (CV) to make the normalized CV plot more comparable.")
    parser.add_argument(
        "--border_penalty",
        action='store_true',
        help="Activate to add a penalty for non-centered peaks.")
    parser.add_argument(
        "--scale_max",
        metavar='float',
        help="Provide a maximum value for the coefficient of variation (CV) plot.")
    parser.add_argument(
        "--seed",
        metavar='int',
        help="Set seed for the optimization scheme.")

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
    ##   MAIN   ##
    ##############

    # Get the outfile name from the input read file.
    outfilename = args.input_bam.split("/")
    outfilename = outfilename[len(outfilename)-1]
    outfilename = outfilename.replace(".bam", "").replace(".bed", "")

    extended_peak_file_name = "{}/peaks_extended.bed".format(args.output_folder)

    # Find maximal peak length and get the number of peaks from the peak file.
    peaks_file = open(args.input_bed, "r")
    max_peak_len = 0
    num_peaks = 0
    for line in peaks_file:
        num_peaks += 1
        data = line.strip("\n").split("\t")
        start = data[1]
        end = data[2]
        length = int(end) - int(start)

        if (length < 25):
            sys.exit("[ERROR] Peak Length has to be at least 25 bases.")

        if (length > max_peak_len):
            max_peak_len = length
    peaks_file.close()

    print("[NOTE] Maximal Peak Length {}.".format(max_peak_len))

    if ( args.length_norm_value ):
        max_peak_len = int(args.length_norm_value)

    if ( max_peak_len <  25 ):
        sys.exit("[ERROR] Maximal Peak Length has to be at least 25 bases.")

    # Extend the peaks to the maximal length if the parameter is set to true.
    if ( args.length_norm ):

        print("[NOTE] Activate length normalization.")

        # Read in chromosome sizes
        chr_sizes_dict = dict()
        chr_sizes_file = open(args.chr_file, "r")
        for line in chr_sizes_file:
            data = line.strip("\n").split("\t")
            if ( data[0] not in chr_sizes_dict ):
                chr_sizes_dict[data[0]] = int(data[1])
        chr_sizes_file.close()

        # Define new coorindate for peaks. Extend to maximal length.
        peaks_file = open(args.input_bed, "r")
        extended_peak_file = open(extended_peak_file_name, "w")

        for line in peaks_file:
            data = line.strip("\n").split("\t")
            start = int(data[1])
            end = int(data[2])
            peak_length = end - start
            extention_left = numpy.round((max_peak_len-peak_length)/2)
            extentions_right = numpy.round((max_peak_len-peak_length)/2)

            # Check if extention left and right make up the max_peak_length, if not,
            # then add or substract randomly either to left or right some extra bases. This happends
            # because of the rounding.
            current_peak_length = extention_left + extentions_right + peak_length
            if ( current_peak_length < max_peak_len ):
                # Set seed if seed is provided.
                if ( args.seed ):
                    numpy.random.seed(int(args.seed))

                if ( numpy.random.randint(low=2, size=1) == 0 ):
                    extention_left +=  max_peak_len - current_peak_length
                else:
                    extentions_right += max_peak_len - current_peak_length

            if ( current_peak_length > max_peak_len):
                # Set seed if seed is provided.
                if (args.seed):
                    numpy.random.seed(int(args.seed))

                if (numpy.random.randint(low=2, size=1) == 0):
                    extention_left -= current_peak_length - max_peak_len
                else:
                    extentions_right -= current_peak_length - max_peak_len

            # Check if extension goes beyond the borders of the chromosome.
            beyond_left = "false"
            if ( (start - extention_left) < 0 ):
                beyond_left = "true"
            beyond_right = "false"
            if ((end + extentions_right) > chr_sizes_dict[data[0]]):
                beyond_right = "true"

            if ( beyond_left == "true" and beyond_right == "false" ):
                extentions_right += extention_left-start
                extention_left = start

            if ( beyond_left == "false" and beyond_right == "true" ):
                extention_left += (end + extentions_right) - chr_sizes_dict[data[0]]
                extentions_right = chr_sizes_dict[data[0]] - end

            if ( beyond_left == "true" and beyond_right == "true" ):
                extention_left = start
                extentions_right = chr_sizes_dict[data[0]] - end

            start = start - extention_left
            end = end + extentions_right

            # A last checkup if peak length is maximum length.
            if ( (end - start) != max_peak_len and not (beyond_left == "true" and beyond_left == "true") ):
                print("[ERROR] Max length of peak not reached.")
                print(data)
                print(start)
                print(end)
                print(end - start)
                print(max_peak_len)

            # Write extended peak to file.
            extended_peak_file.write("{}\t{}\t{}\t{}\n".format(data[0], int(start), int(end), "\t".join(data[3:])))

        peaks_file.close()
        extended_peak_file.close()

    else:
        extended_peak_file_name = args.input_bed

    # Generate Coverage file with bedtools
    coverage_file_name = "{}/{}_coverage.tsv".format(args.output_folder, outfilename)
    sb.Popen("bedtools coverage -a {} -b {} -d -s > {}".format(extended_peak_file_name, args.input_bam,
                                                               coverage_file_name), shell=True).wait()

    print("[NOTE] {} peaks will be evaluated.".format(num_peaks))

    # Dictionaries for the algorithm.
    size_r_peaks_dict = dict()
    prob_success_peaks_dict = dict()
    variance_coverage_peaks_dict = dict()
    num_bp_peaks_dict = dict()
    coordinates_dict = dict()
    strand_dict = dict()
    center_border_diff_left_dict = dict()
    center_border_diff_right_dict = dict()

    # Get the number of lines of the coverage file of bedtools.
    coverage_file = open(coverage_file_name, "r")
    num_coverage_lines = get_line_count(coverage_file)
    coverage_file.close()

    # Calcualte mean and variance of peak coverage profiles
    peak_cov_list = []

    coverage_file = open(coverage_file_name, "r")

    peak_counter = -1
    line_count = 0

    # Go over each line of the bedtools coverage file.
    for line in coverage_file:
        line_count += 1
        data = line.strip("\n").split("\t")
        bp = int(data[len(data)-2])     # bp of the peak
        cov = int(data[len(data)-1])    # Coverage at that bp
        strand = data[5]

        # If the bp == 1 do the negative binomial estimation an start a new peak entry.
        if(bp == 1):
            if( peak_counter != -1 ):
                peak_center = int(num_bp_peaks_dict[peak_counter] / 2)
                peak_center_ext = int(num_bp_peaks_dict[peak_counter] * 0.1)
                # The fit is the alternative version of the NB. But I get the number of successes (r) and the
                # probability of success (p). At least one value needs to be greater than zero, else the estimation
                # makes no sense.
                # The second condition tests for a centered peak. This filters out peaks that have spkied read coverages
                # at the peak ends.
                center_region = peak_cov_list[(peak_center-peak_center_ext-1):(peak_center+peak_center_ext)]
                border_left = peak_cov_list[0:peak_center_ext]
                border_right = peak_cov_list[-peak_center_ext:]
                if ( not all(v == 0 for v in peak_cov_list) and not all(v <= 10 for v in center_region) ):
                    nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                    size_r_peaks_dict[peak_counter] = nb_fit["size"]
                    prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                    variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1-nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
                    center_border_diff_left_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_left)
                    center_border_diff_right_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_right)
                else:
                    size_r_peaks_dict[peak_counter] = 0.0
                    prob_success_peaks_dict[peak_counter] = 0.0
                    variance_coverage_peaks_dict[peak_counter] = 0.0
                    center_border_diff_left_dict[peak_counter] = 0.0
                    center_border_diff_right_dict[peak_counter] = 0.0
            peak_cov_list = []
            peak_cov_list.append(cov)
            peak_counter += 1
            coordinates_dict[peak_counter] = [data[0], data[1], data[2]]
        else:
            peak_cov_list.append(cov)
            num_bp_peaks_dict[peak_counter] = bp
            strand_dict[peak_counter] = strand

            # This condition takes the last line of the coverage file into account. Else I will miss the last entry.
            if ( line_count == num_coverage_lines ):
                peak_center = int(num_bp_peaks_dict[peak_counter] / 2)
                peak_center_ext = int(num_bp_peaks_dict[peak_counter] * 0.1)
                center_region = peak_cov_list[(peak_center - peak_center_ext - 1):(peak_center + peak_center_ext)]
                border_left = peak_cov_list[0:peak_center_ext]
                border_right = peak_cov_list[-peak_center_ext:]
                if (not all(v == 0 for v in peak_cov_list) and not all(v <= 10 for v in center_region)):
                    nb_fit = fnb.fit_nbinom(numpy.array(peak_cov_list))
                    size_r_peaks_dict[peak_counter] = nb_fit["size"]
                    prob_success_peaks_dict[peak_counter] = nb_fit["prob"]
                    variance_coverage_peaks_dict[peak_counter] = (nb_fit["size"] * (1 - nb_fit["prob"])) / (nb_fit["prob"] * nb_fit["prob"])
                    center_border_diff_left_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_left)
                    center_border_diff_right_dict[peak_counter] = numpy.max(center_region) - numpy.max(border_right)
                else:
                    size_r_peaks_dict[peak_counter] = 0.0
                    prob_success_peaks_dict[peak_counter] = 0.0
                    variance_coverage_peaks_dict[peak_counter] = 0.0
                    center_border_diff_left_dict[peak_counter] = 0.0
                    center_border_diff_right_dict[peak_counter] = 0.0

    coverage_file.close()

    filtered_num_peaks = 0
    varcoeff_coverage_peaks_dict = dict()

    if (args.border_penalty):
        print("[NOTE] Activate border penalty.")

    # Calcualte Variantioncoefficient of peak coverage profile.
    for i in range(0, num_peaks):

        # The mean coverage has to be greater than zero or else the VC is not defined.
        if (size_r_peaks_dict[i] > 0):
            # varcoeff_coverage_peaks_dict[i] = 1 / math.sqrt(size_r_peaks_dict[i] * (1 - prob_success_peaks_dict[i]))
            varcoeff_coverage_peaks_dict[i] = math.sqrt((1-prob_success_peaks_dict[i])/size_r_peaks_dict[i])

            if ( args.border_penalty ):
                w1 = 1
                w2 = 1
                if ( center_border_diff_left_dict[i] < 0 ):
                    w1 = 1/abs(center_border_diff_left_dict[i])
                if ( center_border_diff_right_dict[i] < 0 ):
                    w2 = 1/abs(center_border_diff_right_dict[i])
                varcoeff_coverage_peaks_dict[i] = varcoeff_coverage_peaks_dict[i] * w1 * w2

            # Just a safety condition.
            if ( math.isnan(varcoeff_coverage_peaks_dict[i]) ):
                print(varcoeff_coverage_peaks_dict[i])

            # Because the expected number of successes (mean) were estimated, I have to correct the VC based on the
            # changed number of the degrees of freedom.
            # varcoeff_coverage_peaks_dict[i] = varcoeff_coverage_peaks_dict[i] / math.sqrt(num_bp_peaks_dict[i]-1)

            filtered_num_peaks += 1
        else:
            varcoeff_coverage_peaks_dict[i] = -0.01

    print("[NOTE] {} peaks are covered.".format(filtered_num_peaks))

    # Filter our peaks that are completly uncovered.
    filtered_varcoeff_coverage_peaks = []
    for i in varcoeff_coverage_peaks_dict.values():
        if (i >= 0):
            filtered_varcoeff_coverage_peaks.append(i)

    print("[NOTE] Generate CV Plot")

    scale_max = 1.0
    if ( args.scale_max ):
        if ( float(args.scale_max) < 0.0 ):
            sys.exit("[ERROR] Wrong value for scale_max!")
        scale_max = float(args.scale_max)

    # Make vase plot of variationkoefficients.
    f = plt.figure()
    plt.violinplot(filtered_varcoeff_coverage_peaks)
    plt.ylim(0.0, scale_max)
    plt.ylabel('Variationcoefficient of the Peak Profiles')
    f.savefig(args.output_folder + "/VC_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')

    # Normalize all VC so that the scale goes from 0 to 1 which makes a comparison between
    # different profile evaluations easier. Unity-based normalization.
    one = numpy.max(filtered_varcoeff_coverage_peaks)
    zero = numpy.min(filtered_varcoeff_coverage_peaks)
    if ( args.max_norm_value ):
        one = float(args.max_norm_value)

    for i in range(0, len(filtered_varcoeff_coverage_peaks)):
        filtered_varcoeff_coverage_peaks[i] = (filtered_varcoeff_coverage_peaks[i]-zero)/(one-zero)

    print("[NOTE] Generate Normalized CV Plot")

    # Make vase plot of variationkoefficients.
    f = plt.figure()
    plt.violinplot(filtered_varcoeff_coverage_peaks)
    plt.ylim(0.0, 1.0)
    plt.ylabel('Normalized Variationcoefficient of the Peak Profiles')
    f.savefig(args.output_folder + "/Norm_VC_Distribution_{}.pdf".format(outfilename), bbox_inches='tight')

    # Generate the output tabular file.
    print("[NOTE] Generate Output Tabular")
    out_tab_file_name = args.output_folder + "/VC_tab_{}.bed".format(outfilename)
    out_tab_file = open(out_tab_file_name, "w")

    index_sort = numpy.argsort(list(varcoeff_coverage_peaks_dict.values()))[::-1]
    keys_list = list(varcoeff_coverage_peaks_dict.keys())

    for i in index_sort:
        k = keys_list[i]
        coords = coordinates_dict[k]
        # "Chr Start End ID VC Strand bp r p Max_Norm_VC"
        out_tab_file.write("{}\t{}\t{}\tpeak_{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(coords[0], coords[1], coords[2],
                                                                 i+1, varcoeff_coverage_peaks_dict[k], strand_dict[k],
                                                                 num_bp_peaks_dict[k], size_r_peaks_dict[k],
                                                                 prob_success_peaks_dict[k],
                                                                 (varcoeff_coverage_peaks_dict[k]-zero)/(one-zero),
                                                                 center_border_diff_left_dict[k],
                                                                 center_border_diff_right_dict[k]))
    out_tab_file.close()

    # peak_584

    # Sort the tabular file.
    #sb.Popen("sort -r -k5,5 -g {} > {}".format(out_tab_file_name, "{}/VC_tab_{}.bed".format(args.output_folder, outfilename)), shell=True).wait()
    #sb.Popen("rm {}".format(out_tab_file_name), shell=True).wait()

if __name__ == '__main__':
    main()