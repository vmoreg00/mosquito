#!/usr/bin/python
# ls -1 *.maplog > archivos.txt
import sys

if len(sys.argv) < 1:
        print 'Input list of log files is missing.'
        sys.exit(2)
archivo    = sys.argv[1]
names      = []
one 	   = []
overall	   = []

out        = open('summary_map.txt', 'w')
input      = open(archivo, 'r')

for file in input:
	file = file.strip()
	names.append(file)
	log = open(file, 'r')
	for line in log:
		line = line.strip()
		campo = line.split(' ')
		if  len(campo) > 4 and campo[4] == '1':
			one.append(campo[1])
		if  len(campo) > 3 and campo[1] == 'overall':
			overall.append(campo[0])
	log.close()

out.write('Sample\t\tonetime\t\toverall\n')
for i in range(len(names)):
	out.write(names[i] + '\t' + one[i] + '\t' + overall[i] + '\n')
out.close()
