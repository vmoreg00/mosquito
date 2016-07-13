#/bin/bash
#
#                               06/07/16
#                               ---------
#Para cada elemento de la lista
#si no esta en DATADIR
#lo anyadimos al directorio
#A partir  del programa SABRE identificare los individuos segun su barcode
#Debo crear un directorio con barcode
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

if [ ! -d barcodes ]; then mkdir barcodes; fi
#si dentro del  directorio no esta pipiens_barcodes.txt
#muestro cada  linea del documento
#anyado dos campos  los cuales sean uno para los reads uno y  otro para los reads dos de C.pipiens
if [ ! -e barcodes/pipiens_barcodes.txt ]; then
        echo -e 'GGTCGTAAATG\tPipFe1_R1.fastq\tPipFe1_R2.fastq' > barcodes/pipiens_barcodes.txt
        echo -e 'CTAGTACCTG\tPipFe2_R1.fastq\tPipFe2_R2.fastq' >> barcodes/pipiens_barcodes.txt       
        echo -e 'AACTACGGG\tPipFe3_R1.fastq\tPipFe3_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'TCGACGTT\tPipFe6_R1.fastq\tPipFe6_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'GTCAGAGTATG\tPipMa4_R1.fastq\tPipMa4_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'ACTGAGACTG\tPipFe4_R1.fastq\tPipFe4_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'CAGCTCTAG\tPipMa3_R1.fastq\tPipMa3_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'TGATCTCG\tPipMa1_R1.fastq\tPipMa1_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'AGCCATGAATG\tPipMa2_R1.fastq\tPipMa2_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'CCAATGCTTG\tPipMa5_R1.fastq\tPipMa5_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'GATTGCAGG\tPipMa6_R1.fastq\tPipMa6_R2.fastq'>> barcodes/pipiens_barcodes.txt 
        echo -e 'TTGGCATC\tPipFe5_R1.fastq\tPipFe5_R2.fastq'>> barcodes/pipiens_barcodes.txt 
fi
#Ejecutaremos el programa SABRE
#llamamos a sabre
#queremos  como MAX 1 mismatch
#los file tienen barcodes
if [ ! -e PipFe1_R1.fastq ]; then
	sabre pe -m 1 -c \
        	 -f $DATADIR/Cpipiens_S1_L001_R1_001.fastq.gz  \
         	 -r $DATADIR/Cpipiens_S1_L001_R2_001.fastq.gz \
         	 -b barcodes/pipiens_barcodes.txt \
        	 -u $DATADIR/unknown_Cpipiens_R1.fastq.gz \
        	 -w $DATADIR/unknown_Cpipiens_R2.fastq.gz 
fi
LISTA1=(PipFe1 PipFe2 PipFe3 PipFe6 PipMa4 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5)
#Creamos un directorio llamado merged
#lo recorremos y ejecutamos PEAR
#obtendremos un out para casa sample especifico
if [ ! -d merged ]; then mkdir merged; fi
        for i in  0 1 2 3 4 5 6 7 8 9 10 11;do
		if [ ! -e merged/${LISTA1[$i]}'_assembled.fastq ]; then
                	pear  -f ${LISTA1[$i]}'_R1.fastq' \
                      	      -r ${LISTA1[$i]}'_R2.fastq'  \
                      	      -o merged/${LISTA1[$i]} \
                      	      -v 10 \
                      	      -q 15 \
                      	      -j 1 \
                      	      --memory 2G &
		#When we merged the reads, the ones that couldnt be assembled
                # and also were reverse, changed their orientations
                #We use vsearch to do the reverse complementary to get the original one(reverse) 
                #First of all, we do a list with the both forms 
			vsearch --fastx_revcomp merged/${LISTA1[$i]}'.unassembled.reverse.fastq' \
			        --fastqout merged/${LISTA1[$i]}'.unassembled.reverse.reversed.fastq'
		fi
        done
#mismo proceso para molestus
if [ ! -e barcodes/molestus_barcodes.txt ]; then
        echo -e 'GGTCGTAAATG\tmol01_R1.fastq\tmol02_R2.fastq' > barcodes/molestus_barcodes.txt
        echo -e 'GTCAGAGTATG\tMol02_R1.fastq\tMol02_R2.fastq'>> barcodes/molestus_barcodes.txt 
        echo -e 'CTAGTACCTG\tMol03_R1.fastq\tMol03_R2.fastq' >> barcodes/molestus_barcodes.txt
        echo -e 'AACTACGGG\tMol04_R1.fastq\tMol04_R2.fasq'>> barcodes/molestus_barcodes.txt 
        echo -e 'TCGACGTT\tFe6_R1.fastq\tFe6_R2.fastq'>> barcodes/molestus_barcodes.txt 
fi
if [ ! -e Mol01_R1.fastq ]; then
	sabre pe -m 1 -c \
        	 -f $DATADIR/Mol1-5_S1_L001_R1_001.fastq.gz  \
        	 -r $DATADIR/Mol1-5_S1_L001_R2_001.fastq.gz \
        	 -b barcodes/molestus_barcodes.txt \
         	 -u $DATADIR/unknown_Molestus_R1.fastq.gz \
         	 -w $DATADIR/unknown_Molestus_R2.fastq.gz 

LISTA2=(Mol01 Mol02 Mol03 Mol04 Mol05)
#en carpeta merged
if [ ! -d merged ]; then mkdir merged; fi
     for i in  0 1 ;do
	 if [ ! -e merged/${LISTA2[$i]}'.assembled.fastq' ]; then
                pear  -f ${LISTA2[$i]}'_R1.fastq' \
                      -r ${LISTA2[$i]}'_R2.fastq' \
		      -o merged/${LISTA2[$i]} \
                      -v 10 \
                      -q 15 \
                      -j 1 \
                      --memory 2G 
		#When we merged the reads, the ones that couldnt be assembled
		# and also were reverse, changed their orientations
		#We use vsearch to do the reverse complementary to get the original one(reverse) 
		#First of all, we do a list with the both forms	
		 vsearch --fastx_revcomp merged/${LISTA2[$i]}'.unassembled.reverse.fastq' \
                         --fastqout merged/${LISTA2[$i]}'.unassembled.reverse.reversed.fastq'
	fi    
     done
