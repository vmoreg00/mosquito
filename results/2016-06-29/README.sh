#!/bin/bash
#
#				2016-06-29
#				----------
#

DATADIR=`pwd | sed 's/results/data/'`

if [ ! -d $DATADIR ]; then
   mkdir $DATADIR
fi

for i in Mol1-5_S1_L001_R1_001.fastq.gz \
         Mol1-5_S1_L001_R2_001.fastq.gz \
         Cpipiens_S1_L001_R1_001.fastq.gz \
         Cpipiens_S1_L001_R2_001.fastq.gz; do
   if [ ! -e $DATADIR/$i ]; then
      ln -s /data/joiglu/mosquito/$i $DATADIR/$i
   fi
done
