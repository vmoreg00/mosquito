#!/usr/bin/env python3
# correlation test and graph

from matplotlib import pyplot as plt

in_file = "both_qual.sorted.csv"
out = "both_qual.csv"

in_op = open(in_file, "r")
out = open(out, "w+")

in_qual = in_op.readlines()
out.write("#chr\tpos\tfb_alt\tfb_qual\tgatk_alt\tgatk_qual\n")

i = 1
l1 = in_qual[0].strip().split()
fb_qual = []
gatk_qual = []
while i < len(in_qual):
    l2 = in_qual[i].strip().split()
    if l1[2] == l2[2]:
        out.write(l1[1] + "\t" + l1[2] + "\t" + l1[3] + "\t" + l1[4] + "\t"
                  + l2[3] + "\t" + l2[4] + "\n")
        i += 2
        l1 = in_qual[i - 1].strip().split()
    else:
        i += 1
        l1 = l2

in_op.close()
out.close()
