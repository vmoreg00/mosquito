#!/usr/bin/python
# ls -1 *_assembled.log > archivos.txt
import sys
if len(sys.argv) < 1:
	print 'Input assembled log is missing'
	sys.exit(2)

archivo	   = sys.argv[1]
names	   =[]
total	   =[]	 
reads	   =[]
out	   =open('summary_assembled.txt','w')
input	   =open(archivo,'r')

for file in input:
	file= file.strip()
	names.append(file)
	asslog = open(file,'r')
	n=0
	for line in asslog:
		line=line.strip()
		n+=1
		if n==8:
			total.append(line)
		if n==9:
			reads.append(line)
	asslog.close()

out.write('Sample' + '\n' + 'Total reads' +'\n' +'Reads with adapters' + '\n'+ '\n')
for i in range(len(names)):
	out.write(names[i] + '\n' + total[i] + '\n' + reads[1] + '\n'+ '\n')
out.close()
