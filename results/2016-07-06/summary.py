#New file with the log information
#ls -1 *_pear.log > archivos.txt
#calling
import sys
if len(sys.argv)< 1:
        print' Falta output'
        sys.exit()

archivo= sys.argv[1]
#generamos listas utilizadas a continuacion
names=[]
assembled=[]
total=[]
percent=[]
#archivo de salida
out=open('output.txt','w')
#nos recorremos el archivo con todos los file
#para cada file de el archivo con todos los nombre
input= open(archivo,'r')
listasfile= input.readlines()
for file in listasfile:
	file=file.strip()
	names.append(file)
	#para cada linea en open file
	for line in open(file,'r'):
		line=line.strip()
		campo=line.split(' ')

		if len(campo)> 4:
			if campo[0]== 'Assembled' and campo[1]== 'reads':
                        	assembled.append(campo[3])
                       		total.append(campo[5])
                     		percent.append(campo[6])			
#escribimos output
out.write('Sample'+'\t'+'assembled'+'\t'+'total'+'\t'+'percent'+'\n')
for i in range(len(names)-1):
        out.write(names[i] +'\t'+ assembled[i]+ '\t'+total[i]+'\t'+ percent[i]+'\n')

