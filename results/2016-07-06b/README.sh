#!/bin/bash
#
#				2016-07-06b
#				-----------
#
# Some reads may contain adapter sequence, although we do not expect the
# assembled pairs to have any. In any case, it is worth checking and removing
# from downstream analyses those sequences.
MERGED='../2016-07-06/merged'
LIST=(PipFe1 PipFe2 PipFe3 PipFe6 PipMa4 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5 Mol01 Mol02 Mol03 Mol04 Mol05)
PROC=`grep -q -P '^processor' /proc/cpuinfo | wc -l`
if [ $PROC -gt 8 ]; then
   for i in `seq 0 16`; do
      if [ ! -e ${LIST[$i]}_trimmed.fastq ]; then
         # We need to revise this:
         cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -g CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
                  -o ${LIST[$i]}_trimmed.fastq \
                  $MERGED/${LIST[$i]}'.assembled.fastq' > ${LIST[$i]}_assembled.log &
      fi
      if [ ! -e ${LIST[$i]}_R1_trimmed.fastq ]; then
         cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -A CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
                  -o ${LIST[$i]}_R1_trimmed.fastq \
                  -p ${LIST[$i]}_R2_trimmed.fastq \
                  $MERGED/${LIST[$i]}'.unassebled.forward.fastq' \
                  $MERGED/${LIST[$i]}'.unassebled.reverse.fastq' > ${LIST[$i]}.log &
      fi
   done
   wait
else

fi
