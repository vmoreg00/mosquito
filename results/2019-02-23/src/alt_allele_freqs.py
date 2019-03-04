#!/usr/bin/env python3

# Read VCF and extract the alternate and reference alleles frequencies

import sys

if len(sys.argv) != 4:
    print("ERROR: 3 arguments must be provided")
    print("Usage:")
    print("alt_allele_freqs.py <input.VCF> <out_freqs.csv> <out_chr_lengths.csv>")
    sys.exit(1)

# input file
input_vcf = sys.argv[1]
vcf_opened = open(input_vcf, "r")

# output files
freqs_out = sys.argv[2]
lengths = sys.argv[3]
freqs_opened = open(freqs_out, "w+")
lengths_opened = open(lengths, "w+")

# Computations of the alternate allele frequencies
for line in vcf_opened:
    # Extract the header information
    if line[0:2] == "#C":
        populations = line.strip().split("\t")
        header = ["chrom", "pos"] + populations[9:]
        freqs_opened.write("\t".join(header) + "\n")
    # Computes the Alternate allele frequencies
    elif line[0] != "#":
        info = line.strip().split("\t")
        chr = info[0]
        pos = info[1]
        geno = info[9:]
        freqs = []
        for pop in geno:
            geno_info = pop.split(":")
            ref_obs = geno_info[3]
            alt_obs = geno_info[5].split(",")[0]
            if(ref_obs == "." or alt_obs == "."):
                freqs.append("NA")
            elif(alt_obs == "0"):
                freqs.append("0")
            elif(ref_obs == "0"):
                freqs.append("1")
            else:
                freqs.append(str(int(alt_obs) / (int(ref_obs) + int(alt_obs))))
        output = [chr, pos] + freqs
        freqs_opened.write("\t".join(output) + "\n")
    # Extract contig lengths
    elif line[0:8] == "##contig":
        contig = line.strip().split(",")
        name = str(contig[0].split("=")[2])
        length = str(contig[1].split("=")[1].rstrip(">"))
        lengths_opened.write(name + "\t" + length + "\n")

vcf_opened.close()
freqs_opened.close()
lengths_opened.close()

