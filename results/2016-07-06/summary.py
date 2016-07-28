#!/usr/bin/python
#
# This script reads a file with a list of files, which must be
# the standard output of pear, run on several samples. Then, it
# takes some information from each pear output, and writes it down
# in a summary_pear.txt file.
#
# You can create the list of pear output files like this, if they
# have been named correctly:
#
#      ls -1 *_pear.log > archivos.txt
#
# Usage: summary.py <input>

import sys

if len(sys.argv) < 1:
        print 'Input list of pear log files is missing.'
        sys.exit(2)

archivo    = sys.argv[1]
names      = []
assembled  = []
total      = []
percent    = []
out        = open('summary_pear.txt', 'w')
input      = open(archivo, 'r')
listasfile = input.readlines()
input.close()

for file in listasfile:
	file = file.strip()
	names.append(file)
        pearlog = open(file, 'r')
	for line in pearlog:
		line = line.strip()
		campo = line.split(' ')
		if len(campo) > 4 and campo[0] == 'Assembled' and campo[1] == 'reads':
			assembled.append(campo[3])
			total.append(campo[5])
			percent.append(campo[6])
	pearlog.close()

out.write('Sample\tassembled\ttotal\tpercent\n')
for i in range(len(names)):
        out.write(names[i] + '\t' + assembled[i] + '\t'+total[i] + '\t' + percent[i] + '\n')

out.close()

