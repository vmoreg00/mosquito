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
	for line in asslog:
		line=line.strip()
		campo= line.split(' ')
		if campo[0] == 'Total' and campo[1]== 'reads':
			total.append(line)
		if campo[0] == 'Reads' and campo[1] == 'with':
			reads.append(line)
	asslog.close()
out.write('Sample' + '\n' + 'Total reads' +'\n' +'Reads with adapters' + '\n'+ '\n')
for i in range(len(names)):
	out.write(names[i] + '\n' + total[i] + '\n' + reads[i] + '\n'+ '\n')
out.close()
