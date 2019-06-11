#!/usr/bin/env python3

###############################################################################
#                                                                             #
# As PoPoolation2 scripts do not accept VCF files with more than one          #
# popoolation, a script have been developed in order to transform the VCF     #
# mpileup format to synchronize format (requierd by PoPoolation2).            #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2019-06-10                                                            #
#                                                                             #
###############################################################################

import sys

# Argumen checking and parsing
if len(sys.argv) != 3:
    print("ERROR: incorrect arguments")
    print("USAGE:")
    print("./multivcf2sync.py <input.vcf> <output.sync>")
    sys.exit(1)
input = sys.argv[1]
output = sys.argv[2]

# Open files
vcf = open(input, "r")
sync_out = open(output, "w+")

# Convert to synchronize format
for line in vcf:
    # Skip header
    if line[0] == "#":
        continue
    # Extract record information
    line_read = line.strip().split("\t")
    info = {"chr": line_read[0],
            "pos": line_read[1],
            "ref": line_read[3],
            "alt": line_read[4].split(","),
            "pops": [i.split(":")[0] for i in line_read[9:]]
           }
    # Sync format list
    sync = [{"A":0, "C":0, "G":0, "T":0, "N":0, "indel":0}
            for i in range(len(info["pops"]))]
    ## Indel case ...
    if len(info["ref"]) > 1:
        for i in range(len(info["pops"])):
            if info["pops"][i] != ".":
                haplotypes = info["pops"][i].split("/")
                sync[i]["indel"] = len(haplotypes(
    ## ... SNP case
    else:
        for i in range(len(info["pops"])):
            haplotypes = info["pops"][i].split("/")
            sync[i][info["ref"]] += haplotypes.count("0")
            hapl_number = 1
            for j in info["alt"]:
                sync[i][j] += haplotypes.count(str(hapl_number))
                hapl_number += 1
    # Write output
    sync_info = [info["chr"], info["pos"], info["ref"]]
    for pop in sync:
        alt_count = [str(i) for i in list(pop.values())]
        alt_count = ":".join(alt_count)
        sync_info.append(alt_count)
    sync_out.write("\t".join(sync_info) + "\n")

# Close files
vcf.close()
sync_out.close()
