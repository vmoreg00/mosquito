#!/bin/bash
#
#				06-07-16b
#				----------
#A partir de las reads ensambladas
#quitaremos el adaptador con el programa cutadapt
#nos la recorremos
LISTA1=( PipFe1 Mol01 PipFe2 Mol02 PipFe3 Mol03 PipFe6 Mol04 PipMa4 Mol05 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5)
#Para cada elemento de la lista
#ejecutaremos cut adapt otorgandole los dos adaptadores
#primero lo ejecutaremos para las reads que fueron ensambladas
#posteriormente a las read 1 y 2 de las no ensambladas, donde se espera tener mas
for i in  0 1 ;do	
		cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                 	 -g CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
                	 -o ${LISTA1[$i]}_trimmed.fastq \
                 	 /home/student/Documents/mosquito/results/06-07-16/merged/${LISTA1[$i]}'.assebled.fastq' > cutadapt.log &

		cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
			 -A CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
			 -o ${LISTA1[$i]}_R1_trimmed.fastq \
			 -p ${LISTA1[$i]}_R2_trimmed.fastq \
			 /home/student/Documents/mosquito/results/06-07-16/merged/${LISTA1[$i]}'.unassebled.forward.fastq' \
			 /home/student/Documents/mosquito/results/06-07-16/merged/${LISTA1[$i]}'.unassebled.reverse.fastq' > cutadapt.log &
done
wait
