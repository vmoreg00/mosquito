#/bin/bash
#
#                               06/07/16
#                               ---------
#Para cada elemento de la lista
#si no esta en DATADIR
#lo anyadimos al directorio
#A partir  del programa SABRE identificare los individuos segun su barcode
#Debo crear un directorio con barcode
if [ ! -e /home/student/Documents/barcodes/pipiens_barcodes.txt ]; then
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
sabre pe -m 1 -c \
         -f /home/student/Documents/data/pipiens_R1.fastq \
         -r /home/student/Documents/data/pipiens_R2.fastq  \
         -b barcodes/pipiens_barcodes.txt \
         -u /home/student/Documents/data/unknown_Cpipiens_R1.fastq.gz \
         -w /home/student/Documents/data/unknown_Cpipiens_R2.fastq.gz 

LISTA1=(Fe1 Fe2 Fe3 Fe6 Ma4 Fe4 Ma3 Ma1 Ma2 Ma5 Ma6 Fe5)
#Creamos un directorio llamado merged
#lo recorremos y ejecutamos PEAR
#obtendremos un out para casa sample especifico
if [ ! -d merged ]; then mkdir merged; fi
     for i in  0 1 2 3 ;do
	#si no esta presente la muestra en merged
	#ejecutaremos pear para que aparezca
                pear  -f /home/student/Documents/mosquito/results/06-07-16/${LISTA1[$i]}'_R1.fastq' \
                      -r /home/student/Documents/mosquito/results/06-07-16/${LISTA1[$i]}'_R2.fastq' \
                      -o merged/${LISTA1[$i]} \
		      -v 10 \
                      -q 15 \
                      -j 1 \
                      --memory 2G 
		#When we merged the reads, the ones that couldnt be assembled
                # and also were reverse, changed their orientations
                #We use vsearch to do the reverse complementary to get the original one(reverse) 
                #First of all, we do a list with the both forms 
		 vsearch --fastx_revcomp merged/${LISTA1[$i]}'.unassembled.reverse.fastq' \
                         --fastqout merged/${LISTA1[$i]}'.unassembled.reverse.reversed.fastq'
     done
#mismo proceso para molestus
if [ ! -e /home/student/Documents/barcodes/molestus_barcodes.txt ]; then
        echo -e 'GGTCGTAAATG\tMol01_R1.fastq\tMol01_R2.fastq' > barcodes/molestus_barcod$
        echo -e 'GTCAGAGTATG\tMol02_R1.fastq\tMol02_R2.fastq'>> barcodes/molestus_barcod$
        echo -e 'CTAGTACCTG\tMol03_R1.fastq\tMol03_R2.fastq' >> barcodes/molestus_barcod$
        echo -e 'AACTACGGG\tMol04_R1.fastq\tMol04_R2.fastq'>> barcodes/molestus_barcodes$
        echo -e 'TCGACGTT\tMol05_R1.fastq\tMol05_R2.fastq'>> barcodes/molestus_barcodes.$
fi
sabre pe -m 1 -c \
         -f /home/student/Documents/data/molestus_R1.fastq \
         -r /home/student/Documents/data/molestus_R2.fastq  \
         -b /home/student/Documents/barcodes/molestus_barcodes.txt \
         -u /home/student/Documents/data/unknown_molestus_R1.fastq.gz \
         -w /home/student/Documents/data/unknown_molestus_R2.fastq.gz

LISTA2=(Mol01 Mol02 Mol03 Mol04 Mol05)
#en carpeta merged
if [ ! -d merged ]; then mkdir merged; fi
     for i in  0 1 ;do
                pear  -f /home/student/Documents/mosquito/results/06-07-16/${LISTA2[$i]}'_R1.fastq' \
                      -r /home/student/Documents/mosquito/results/06-07-16/${LISTA2[$i]}'_R2.fastq' \
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
    done

		
