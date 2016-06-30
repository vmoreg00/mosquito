#!/bin/bash
#
#				30-06-16
#				----------
#
#We want to merge the reads so
#We are gonna to use PEAR
#a fast and accurate Illumina Paired-End reAd mergeR

DATADIR=`pwd | sed 's/results/data/'`
if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi
#para ejecutar el programa creamos un directorio
if [ ! -d $merged]; then mkdir $merged; fi

#barcode for molestus
if [ ! -d barcodes ]; then mkdir barcodes; fi
#si dentro del  directorio no esta pipiens_barcodes.txt
#muestro cada  linea del documento
#anyado dos campos  los cuales sean uno para los reads uno y  otro par
#lo ejecutamos

if [ ! -e barcodes/moletus_barcodes.txt ]; then
	echo -e 'GGTCGTAAATG\tmol01_R1.fastq\tmol02_R2.fastq' > barcodes/molestus_barcodes.tx
	echo -e 'GTCAGAGTATG\tMol02_R1.fastq\tMol02_R2.fastq'>> barcodes/molestus_barcodes.txt 
	echo -e 'CTAGTACCTG\tMol03_R1.fastq\tMol03_R2.fastq' >> barcodes/molestus_barcodes.txt
	echo -e 'AACTACGGG\tMol04_R1.fastq\tMol04_R2.fasq'>> barcodes/molestus_barcodes.txt 
	echo -e 'TCGACGTT\tFe6_R1.fastq\tFe6_R2.fastq'>> barcodes/molestus_barcodes.txt 
fi


