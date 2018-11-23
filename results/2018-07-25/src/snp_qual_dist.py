#!/usr/bin/env python3

###############################################################################
#                                                                             #
# snp_qual_distribution.py                                                    #
#                                                                             #
# Determines the frequency distribution of the quality of the SNPs included   #
# in a vcf file                                                               #
#                                                                             #
# Usage: python3 snp_qual_distribution.py <file.vcf> <output>                 #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2018-10-18                                                            #
#                                                                             #
###############################################################################

import sys
import matplotlib.pyplot as plt

def read_vcf(filename):
    """
    Reads a VCF file and extracts the quality of the SNPs only.

    Args:
        filename (str) Name of the VCF file to be read

    Returns:
        list: Qualities of the SNPs in VCF file
    """

    vcf_opened = open(filename, "r")

    # Extracts the quality of SNPs
    qual = []
    for line in vcf_opened:
        if line[0] != "#":
            reg = line.strip().split("\t")
            alt = reg[4].split(",")
            if len(alt[0]) == 1:
                qual.append(float(reg[5]))
    vcf_opened.close()
    return qual


def write_qual(qual, freq, filename):
    """
    Writes the absolute frequencies of the SNPs in a CSV file

    Args:
        qual (list) List of the qualities of the SNPs
	freq (list) Frequency of the intervals of qual
        filename (str) Name of the CSV file to be wroten

    Returns:
        CSV file
    """

    output = open(filename, "w+")
    for i in range(len(freq)):
        interval = "[" + str(qual[i]) + ", " + str(qual[i+1]) + ")"
        frq = str(freq[i])
        output.write(interval + "\t" + frq + "\n")
    output.close


def threshold(qual, thr):
    """
    Modifies the list of qualities and transform any value larger than the
    threshold into the threshold value.

    Args:
        qual (list) List of the qualities of the SNPs
        thr (int | float) Threshold above which values will be transformed

    Returns:
        Modified list
    """

    for i in range(len(qual)):
        if qual[i] > thr:
            qual[i] = thr


def histogram(qual, nbins, filename):
    """
    Creates the distribution of frequencies of the qualities of the SNPs.

    Args:
        qual (list) List of the qualities of the SNPs
        nbins (int) Number of bins to represent in the histogram
        filename (str) Name of the output file in which the plot will be saved

    Returns:
        Image and CSV file

    """

    name = filename.split("/")[-1]
    name = name.split(".")[0]

    n, bins, patches = plt.hist(qual, nbins)
    plt.xlabel("Quality (phred score)")
    plt.ylabel("Absolute frequency")
    plt.title(name)
    plt.savefig(filename + ".png")

    write_qual(bins, n, filename + ".csv")


def main():
    # Information to the user
    if len(sys.argv) != 3:
        print("\nDetermines the frequency distribution of the quality of the\n"
              + "SNPs included in a vcf file.\n\n"
              + "Usage: snp_qual_distribution.py <file.vcf> <output>")
        sys.exit(1)

    # Argument parsing
    vcf = str(sys.argv[1])
    output = str(sys.argv[2])

    # File parsing
    qual = read_vcf(vcf)

    # Histogram creation
    threshold(qual, 400)
    histogram(qual, 100, output)


if __name__ == "__main__":
    main()
