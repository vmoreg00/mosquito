#!/usr/bin/env python3

###############################################################################
#                                                                             #
# This script computes the mean and standard deviation across all Fst values  #
# obtained in the sliding windows analsyis with PoPoolation.                  #
#                                                                             #
# Only windows with more than 10 SNPs are considered.                         #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2019-06-10                                                            #
#                                                                             #
###############################################################################
import sys
from math import sqrt

# Argument checking and parsing
if len(sys.argv) != 3:
    print("ERROR: incorrect arguments")
    print("USAGE:")
    print("./avg-fst.py <input.vcf> <output_preffix>")
    sys.exit(1)
input = sys.argv[1]
output = sys.argv[2]

# Initialize mean and SD matrices
fst_avg = [[0 for j in range(8)] for i in range(8)]
fst_sd = [[0 for j in range(8)] for i in range(8)]

# Arithmetic mean
## Sum the terms
sync = open(input, "r")
counter = 0
for line in sync:
    info = line.strip().split("\t")
    if int(info[2]) >= 10:
        counter += 1
        for pair in info[5: ]:
            pair_list = pair.split("=")
            pops = pair_list[0].split(":")
            p1 = int(pops[0]) - 1
            p2 = int(pops[1]) - 1
            fst = float(pair_list[1])
            fst_avg[p1][p2] += fst
sync.close()

## Divide by sample size
for i in range(len(fst_avg)):
    for j in range(len(fst_avg[i])):
        fst_avg[i][j] /= counter

## Save the results
avg_out = open(output + "_mean.csv", "w+")
avg_out.write("\t" + "\t".join([str(i + 1) for i in range(8)]) + "\n")
for i in range(len(fst_avg)):
    avg_out.write(str(i + 1) + "\t" \
                  + "\t".join(str(j) for j in fst_avg[i]) + "\n")
avg_out.close()


# Standard deviation
## Sum the terms
sync = open(input, "r")
for line in sync:
    info = line.strip().split("\t")
    if int(info[2]) >= 10:
        for pair in info[5: ]:
            pair_list = pair.split("=")
            pops = pair_list[0].split(":")
            p1 = int(pops[0]) - 1
            p2 = int(pops[1]) - 1
            fst = float(pair_list[1])
            fst_sd[p1][p2] += (fst - fst_avg[p1][p2]) ** 2
sync.close()

## SD formula
for i in range(len(fst_sd)):
    for j in range(len(fst_sd[i])):
        fst_sd[i][j] = sqrt(fst_sd[i][j] / (counter - 1))

## Save the results
sd_out = open(output + "_sd.csv", "w+")
sd_out.write("\t" + "\t".join([str(i + 1) for i in range(8)]) + "\n")
for i in range(len(fst_sd)):
    sd_out.write(str(i + 1) + "\t" \
                 + "\t".join(str(j) for j in fst_sd[i]) + "\n")
sd_out.close()
