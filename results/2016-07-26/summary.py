#!/usr/bin/python
#ls -1 stats/*.txt >archivos.txt
import sys

if len(sys.argv) < 1:
	print 'Cluster file is missing'
	sys.exist(2)

files	= open(sys.argv[1],'r')
out    = open('summary_clust.txt','w')

total	= {}
dptme	= {}
d1tot	= {}
d1me	= {}
lista=[]
for archivo in files:
	archivo=archivo.strip()
	clust=open(archivo,'r')
	n=0
	for line in clust:
		line=line.strip()
		n+=1
		if n>2 and n< 15:
			campo=line.split('\t')
			#Dictionaries: key + value
			total[archivo + '.' + campo[0]] = 	campo[2]
			lista.append(archivo + '.' + campo[0])
			dptme[archivo + '.' + campo[0]] = 	campo[3]
			d1tot[archivo + '.' + campo[0]] =	campo[4]
			d1me[archivo  + '.' + campo[0]] =	campo[5]
	clust.close()
lista.sort()
out.write('Wclust\tmuestra\ttotal\tdptme\td>1.tot\td>1.me\n')
for i in lista:	
	z1=i.split('.')
        wclust= z1[1]
        muestra= z1[3]
	out.write( wclust + '\t' + muestra + '\t' +  total[i] + '\t'  + dptme[i] + '\t'  + d1tot[i] + '\t' + d1me[i] + '\n')

out.close()

