#!/usr/bin/python
# ls -1 *_paired.log > archivos_pe.txt
import sys
if len(sys.argv) < 1:
        print 'Input paired log is missing'
        sys.exit(2)
archivo    = sys.argv[1]
names      =[]
r1     	   =[]
r2         =[]
short      =[]
out        =open('summary_paired.txt','w')
input      =open(archivo,'r')

for file in input:
        file= file.strip()
        names.append(file)
        pairlog = open(file,'r')
	n=0
        for line in pairlog:
		line=line.strip()
		campo= line.split(' ')
		n+=1
                line=line.strip()
		if n> 7 and n<12:
                	if n == 9:
	                        r1.append(line)
        	        if n == 10:
                                r2.append(line)
			if n== 11:
				if campo[0] == 'Pairs' and campo[1]== 'that':
					short.append(line)
				else:
					short.append('------')
	pairlog.close()
out.write('Sample' + '\n' + 'Read1 with adapter' + '\n'+ 'Read 2 with adapter' + '\n' + '*Pairs that were too short'+ '\n'+'\n')
for i in range(len(names)):
        out.write(names[i] + '\n' + r1[i] +'\n' + r2[i] + '\n' + short[i]+ '\n'+ '\n')
out.close()
