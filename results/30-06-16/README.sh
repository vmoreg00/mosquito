
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

#Para ejecutar PEAR
#realizaremos uns lista con todas las muestras
#De pipiens
LISTA=(Fe1, Fe2, Fe3, Fe6, Ma4, Fe4, Ma3, Ma1, Ma2, Ma5, Ma6, Fe5)

#Creamos un directorio llamado merged
#para cada elemento de la lista
#lo recorremos y ejecutamos PEAR
if [ ! -d $merged]; then mkdir $merged; fi
	for i in  0,1,2,3,4,5,6,7,8,9,10,11;do
		pear  -f ${LISTA[$i]}_R1.fastq \
      		      -r ${LISTA[$i]}_R2.fastq  \
                      -o $merged/'out_merged${LISTA[$i]}_Cpipiens.fastq.gz' \
                      -v 10 \
                      -q 15 \
                      -j 1 \
                      --memory 800M &

	done

