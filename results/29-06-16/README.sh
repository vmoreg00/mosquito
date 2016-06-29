#/bin/bash
#
#				29/06/16
#				----------

#We have 2 types of mosquito, C.pipiens and molestus. We have seqs of both
#Each one with 2 different docs which contains reads
#First of all, we are going to evaluate the quality of the sequences
#and the effect of the coberage

#Defino variable DATADIR donde sustituyo por data
#si no existe la creo

DATADIR=`pwd | sed 's/results/data/'` 
if [ ! -d $DATADIR ]; then
	mkdir $DATADIR
fi
#Para cada elemento de la lista
#si no esta en DATADIR
#lo anyadimos al directorio
for i in Mol1-5_S1_L001_R1_001.fastq.gz \
         Mol1-5_S1_L001_R2_001.fastq.gz \
         Cpipiens_S1_L001_R1_001.fastq.gz \
         Cpipiens_S1_L001_R2_001.fastq.gz; do
   if [ ! -e $DATADIR/$i ]; then
      ln -s /data/joiglu/mosquito/$i $DATADIR/$i
   fi
done

#A partir  del programa SABRE identificare los individuos segun su barcode
#Debo crear un directorio con barcode

if [ ! -d barcodes ]; then mkdir barcodes; fi
#si dentro del  directorio no esta pipiens_barcodes.txt
#muestro cada  linea del documento
#anyado dos campos  los cuales sean uno para los reads uno y  otro para los reads dos de C.pipiens
if [ ! -e barcodes/pipiens_barcodes.txt ]; then
	echo -e 'GGTCGTAAATG \t Fe1_R1.fastq\t Fe1_R2.fastq' > barcodes/pipiens_barcodes.txt
	echo -e 'CTAGTACCTG \t Fe2_R1.fastq\t Fe2_R2.fastq' >> barcodes/pipiens_barcodes.txt	
	echo -e 'AACTACGGG \t Fe3_R1.fastq\t Fe3_R2.fasq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'TCGACGTT \t Fe6_R1.fastq\t Fe6_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'GTCAGAGTATG \t Ma4_R1.fastq\t Ma4_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'ACTGAGACTG \t Fe4_R1.fastq\t Fe4_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'CAGCTCTAG \t Ma3_R1.fastq\t Ma3_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'TGATCTCG \t Ma1_R1.fastq\t Ma1_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'AGCCATGAATG \t Ma2_R1.fastq\t Ma2_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'CCAATGCTTG \t Ma5_R1.fastq\t Ma5_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'GATTGCAGG \t Ma6_R1.fastq\t Ma6_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'TTGGCATC \t Fe5_R1.fastq\t Fe5_R2.fastq'>> barcodes/pipiens_barcodes.txt 
fi
