#!/bin/bash
#
#				2016-07-06b
#				-----------
#
# Some reads may contain adapter sequence, although we do not expect the
# assembled pairs to have any. In any case, it is worth checking and removing
# from downstream analyses those sequences.
#
# The adapter sequences are the following:
#
#    Adapter 1: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA
#    Adapter 2: CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA
#
# Adapter 1 is ligated on the 5' end of the first or forward read, so that its
# reverse complementary may appear on the 3' end of the second or reverse read.
# Adapter 2 is ligated on the 5' end of the second read, so that its reverse
# complementary may appear on the 3' end of the first read. See if this helps:
#
#                     read 1
#      Adap. 1  ------------------>       Barcode.
#   5-----------+---+---------------------+---+--------------3
#    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#   3-----------+---+---------------------+---+--------------5
#              Barcode       <-----------------     Adap. 2
#                                   read 2
#
# The reverse complementary sequences are these:
#    RevComp(Adapter 1): TCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#    RevComp(Adapter 2): TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG

MERGED='../2016-07-06/merged'
LIST=(PipFe1 PipFe2 PipFe3 PipFe6 PipMa4 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5 Mol01 Mol02 Mol03 Mol04 Mol05)
PROC=`grep -P '^processor' /proc/cpuinfo | wc -l`
if [ $PROC -gt 34 ]; then
   for i in `seq 0 16`; do
      if [ ! -e ${LIST[$i]}_setrimmed.fastq ]; then
         # We need to revise this:
         cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG \
                  -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -o ${LIST[$i]}_setrimmed.fastq \
                  -m 30 \
                  $MERGED/${LIST[$i]}'.assembled.fastq' > ${LIST[$i]}_assembled.log &
      fi
      if [ ! -e ${LIST[$i]}_R1_trimmed.fastq ]; then
         cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG \
                  -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -A TCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                  -G CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
                  -o ${LIST[$i]}_R1_trimmed.fastq \
                  -p ${LIST[$i]}_R2_trimmed.fastq \
                  -m 30 \
                  $MERGED/${LIST[$i]}'.unassembled.forward.fastq' \
                  $MERGED/${LIST[$i]}'.unassembled.reverse.fastq' > ${LIST[$i]}_paired.log &
      fi
   done
   wait
else
   for j in 0 3 6 9 12; do
      for i in `seq $j $(( j + 2 ))`; do
      if [ ! -e ${LIST[$i]}_setrimmed.fastq ]; then
         # We need to revise this:
         cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG \
                  -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -o ${LIST[$i]}_setrimmed.fastq \
                  -m 30 \
                  $MERGED/${LIST[$i]}'.assembled.fastq' > ${LIST[$i]}_assembled.log &
      fi
      if [ ! -e ${LIST[$i]}_R1_trimmed.fastq ]; then
         cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG \
                  -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -A TCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                  -G CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
                  -o ${LIST[$i]}_R1_trimmed.fastq \
                  -p ${LIST[$i]}_R2_trimmed.fastq \
                  -m 30 \
                  $MERGED/${LIST[$i]}'.unassembled.forward.fastq' \
                  $MERGED/${LIST[$i]}'.unassembled.reverse.fastq' > ${LIST[$i]}_paired.log &
      fi
      done
      wait
   done
   for i in 15 16; do
      if [ ! -e ${LIST[$i]}_setrimmed.fastq ]; then
         # We need to revise this:
         cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG \
                  -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -o ${LIST[$i]}_setrimmed.fastq \
                  -m 30 \
                  $MERGED/${LIST[$i]}'.assembled.fastq' > ${LIST[$i]}_assembled.log &
      fi
      if [ ! -e ${LIST[$i]}_R1_trimmed.fastq ]; then
         cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCTATGTATCTCGTATGCCGTCTTCTGCTTG \
                  -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
                  -A TCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                  -G CAAGCAGAAGACGGCATACGAGATACATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGA \
                  -o ${LIST[$i]}_R1_trimmed.fastq \
                  -p ${LIST[$i]}_R2_trimmed.fastq \
                  -m 30 \
                  $MERGED/${LIST[$i]}'.unassembled.forward.fastq' \
                  $MERGED/${LIST[$i]}'.unassembled.reverse.fastq' > ${LIST[$i]}_paired.log &
      fi
   done
   wait
fi

if [ ! -e summary_assembled.txt ]; then
   if [ ! -e archivos_se.txt ]; then
      ls -1 *_assembled.log > archivos_se.txt
   fi
   python summary.py archivos_se.txt
   rm archivos_se.txt
fi

if [ ! -e summary_paired.txt ]; then
   if [ ! -e archivos_pe.txt ]; then
      ls -1 *_paired.log > archivos_pe.txt
   fi
   python summary2.py archivos_pe.txt
   rm archivos_pe.txt
fi
