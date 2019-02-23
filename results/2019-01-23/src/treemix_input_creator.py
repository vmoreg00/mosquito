#!/usr/bin/env python3

"""
Creates the treemix input file from a modifyied VCF (only the genotypes
without header).

Name: treemix_input_creator.py
Author: Victor Moreno-Gonzalez
Date: 2019-01-23
"""

import sys

# Check the arguments
if len(sys.argv) != 3:
    print("ERROR: 2 arguments required and %d were given" %
          (len(sys.argv) - 1))
    print("USAGE: treemix_input_creator.py <INPUT.FILE> <OUTPUT.FILE>")
    sys.exit(1)

# Global variables
in_file = sys.argv[1]
out_file = sys.argv[2]
header = "0.quin 1.molA1 2.pipA4 3.molM1 4.torM2 5.molS1 6.mixS2 7.mixS3 8.torM4"

# Open the files
in_opened = open(in_file, "r")
out_opened = open(out_file, "w+")

# Parse the information of allele coutns
out_opened.write(header + "\n")
for line in in_opened:
    info = line.strip().split()
    SNP_freqs = ["2,0", ]
    for field in info:
        if field == ".": # Empty field
            SNP_freqs.append("0,0")
        else:            # Non-empty field
            SNP_pop_i = field.strip().split(":")
            ref_count = SNP_pop_i[3]
            alt_count = SNP_pop_i[5]

            if ref_count == ".":
                ref_count = "0"
            elif "," in ref_count:
                ref_count = str(ref_count.split(",")[0])

            if alt_count == ".":
                alt_count = "0"
            elif "," in alt_count:
                alt_count = str(alt_count.split(",")[0])

            if ref_count == "0" and ref_count == "0":
                add = False
                break
            else:
                add = True
                SNP_freqs.append(",".join((ref_count, alt_count)))
    if add:
        out_opened.write(" ".join(SNP_freqs) + "\n")

# Close files
in_opened.close()
out_opened.close()
