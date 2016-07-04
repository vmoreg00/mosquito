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
	echo -e 'GGTCGTAAATG\tFe1_R1.fastq\tFe1_R2.fastq' > barcodes/pipiens_barcodes.txt
	echo -e 'CTAGTACCTG\tFe2_R1.fastq\tFe2_R2.fastq' >> barcodes/pipiens_barcodes.txt	
	echo -e 'AACTACGGG\tFe3_R1.fastq\tFe3_R2.fasq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'TCGACGTT\tFe6_R1.fastq\tFe6_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'GTCAGAGTATG\tMa4_R1.fastq\tMa4_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'ACTGAGACTG\tFe4_R1.fastq\tFe4_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'CAGCTCTAG\tMa3_R1.fastq\tMa3_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'TGATCTCG\tMa1_R1.fastq\tMa1_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'AGCCATGAATG\tMa2_R1.fastq\tMa2_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'CCAATGCTTG\tMa5_R1.fastq\tMa5_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'GATTGCAGG\tMa6_R1.fastq\tMa6_R2.fastq'>> barcodes/pipiens_barcodes.txt 
	echo -e 'TTGGCATC\tFe5_R1.fastq\tFe5_R2.fastq'>> barcodes/pipiens_barcodes.txt 
fi
#Ejecutaremos el programa SABRE
#llamamos a sabre
#queremos  como MAX 1 mismatch
#los file tienen barcodes
sabre pe -m 1 -c \
	 -f $DATADIR/Cpipiens_S1_L001_R1_001.fastq.gz  \
         -r $DATADIR/Cpipiens_S1_L001_R2_001.fastq.gz \
	 -b barcodes/pipiens_barcodes.txt \
	 -u $DATADIR/unknown_Cpipiens_R1.fastq.gz \
	 -w $DATADIR/unknown_Cpipiens_R2.fastq.gz 


