#!/usr/bin/python
# ls -1 *_paired.log > archivos_pe.txt
import sys
if len(sys.argv) < 1:
        print 'Input paired log is missing'
        sys.exit(2)
archivo    = sys.argv[1]
names      =[]
rpairs     =[]
bp         =[]
filtered   =[]
out        =open('summary_paired.txt','w')
input      =open(archivo,'r')

for file in input:
        file= file.strip()
        names.append(file)
        pairlog = open(file,'r')
        n=0
        for line in pairlog:
                line=line.strip()
                n+=1
                if n==8:
                        rpairs.append(line)
                if n==13:
                        bp.append(line)
		if n==16:
			filtered.append(line)
        pairlog.close()

out.write('Sample' + '\n' + 'Total read pairs processed' +'\n' +'Total basepairs processed' + '\n'+ 'Total filtered'+ '\n'+'\n')
for i in range(len(names)):
        out.write(names[i] + '\n' + rpairs[i] + '\n' + bp[i] + '\n'+ filtered[i]+ '\n'+ '\n')
out.close()
