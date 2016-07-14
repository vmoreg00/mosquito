#!/bin/bash
#
#				08-07-16
#				----------
#Check why the mapping wasnt very good
# we are gonnna use Bowtie2
if [ ! -e reference.fa ]; then
	wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz
	mv Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz reference.fa.gz
	gunzip reference.fa.gz
fi
#Creo una lista con los nombres
#SI no existe un directorio de mapeo lo creo
#para cada muestra ejecutare bowie
LISTA=(Fe1 Fe2 Fe3 Fe4 Fe5 Fe6 Ma1 Ma2 Ma3 Ma4 Ma5 Ma6 Mol01 Mol02 Mol03 Mol04 Mol05)
if [ ! -d mappeo ]; then mkdir mappeo; fi
if [ ! -e culex.1.bt2 ]; then
 #Hago referencia indexada y ordenada
                 bowtie2-build reference.fa culex
fi
for i in 0 1 2 3 4 5 6 7 8 9 10; do
#Para cada una de las muestras
#si no existe el archivo final sam de mapeo
#ejecuto bowtie
	if [ ! -e mappeo/${LISTA[$i]}'r1'.sam ]; then
		#Ejecuto
		bowtie2		 --local \
				 --very-sensitive \
				 -x culex \
	 	 	 	 -U /home/student/Documents/mosquito/results/30-06-16/${LISTA[$i]}'_R1.fastq' \
		 	 	 -S mappeo/${LISTA[$i]}'r1'.sam \
				 --rg-id ${LISTA[$i]}
	fi
	if [ ! -e mappeo/${LISTA[$i]}'r2'.sam ]; then
                #Ejecuto
                bowtie2          --local \
                                 --very-sensitive \
                                 -x culex \
                                 -U /home/student/Documents/mosquito/results/30-06-16/${LISTA[$i]}'_R2.fastq' \
                                 -S mappeo/${LISTA[$i]}'r2'.sam \
                                 --rg-id ${LISTA[$i]}
        fi
done
#Una vez obtenidos los archivos en .sam, convertirlos a formato bam (binario de sam)
#usaremos el programa sam tools
#como tenemos header el comando usado en samtools es -bS
#Si no existe directorio BAM lo creamos
#Para cada elemento de la lista si no existe ya el archivo
#convertiremos a bam
if [ ! -d BAM ]; then mkdir BAM; fi
for i in 0 1 2 3 4 5 6 7 8 9 10; do
	if [ ! -e BAM/${LISTA[$i]}'r1'.bam ]; then
		samtools view -bS mappeo/${LISTA[$i]}'r1'.sam > BAM/${LISTA[$i]}'r1'.bam
	fi
	if [ ! -e BAM/${LISTA[$i]}'r2'.bam ]; then
                samtools view -bS mappeo/${LISTA[$i]}'r2'.sam > BAM/${LISTA[$i]}'r2'.bam
	fi
done	

