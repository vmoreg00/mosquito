#!/usr/bin/env python3

###############################################################################
#                                                                             #
# filter_vcf_positions.py                                                     #
#                                                                             #
# Removes the positions of a vcf file that are not included in a position     #
# file. The motivation of this script is because VCFTOOLS cannot handle       #
# poliploydi.                                                                 #
#                                                                             #
# Usage: python3 filter_vcf_positions.py <in.vcf> <pos.csv> <out.vcf>         #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2018-10-28                                                            #
#                                                                             #
###############################################################################

import sys

def filter_vcf(vcf, filter, out):
    """
    Filter the in VCF file according to the positions included in filter.

    Args:
        vcf (str) input VCF file
        filter (str) tab-separated file with two columns: chromosome and
                     position
        out (str) output VCF file
    Returns:
        kept (int) Number of variants kept
        removed (int) Number of variants removed
    """

    # Files opening
    vcf_open = open(vcf, "r")
    pos_open = open(filter, "r")
    out_open = open(out, "w+")

    # Files parsing
    kept = 0
    removed = 0
    pos_info = pos_open.readline()
    pos_info = pos_info.strip().split()
    for line in vcf_open:
        if line[0] == "#":
            out_open.write(line)
        else:
            vcf_info = line.strip().split()
            if vcf_info[0] == pos_info[0] and vcf_info[1] == pos_info[1]:
                out_open.write(line)
                pos_info = pos_open.readline()
                pos_info = pos_info.strip().split()
                kept += 1
            else:
                removed += 1
        if len(pos_info) == 0:
            break

    vcf_open.close()
    pos_open.close()
    out_open.close()

    return kept, removed

def main():
    # Information to the user
    if len(sys.argv) != 4:
        print("\nRemoves the positions of a vcf file that are not included",
              "\nin a position file.\n",
              "\nUsage: python3 filter_vcf_positions.py <in.vcf> <pos> <out.vcf>",
              "\n\tin.vcf:  VCF input file",
              "\n\tpos:     tab-separated file with two columns: chromosome and",
              "\n\t         position to keep",
              "\n\tout.vcf: VCF output file")
        sys.exit(1)

    # Argument parsing
    vcf = str(sys.argv[1])
    pos = str(sys.argv[2])
    out = str(sys.argv[3])

    # Files parsing
    k, r = filter_vcf(vcf, pos, out)

    # Final info
    print("%d variants were kept (%d removed) for file %s" % (k, r, vcf))

if __name__ == "__main__":
    main()
