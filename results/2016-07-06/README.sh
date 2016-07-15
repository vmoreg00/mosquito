#!/bin/bash
#
#                               2016-07-06
#                               ----------
#
# Here, we are running first SABRE to demultiplex the original fastq files,
# and then PEAR to merge the paired reads that may overlap.
#
# In order to run this with a reduced dataset, just to test the scripts, I
# define an alternative path for the data.
#
# ========================================================================
# CUSTOMIZE PATHS BELOW!

#SOURCEDIR='/home/student/Documents/data'
#FILENAME=( pipiens_R1.fastq pipiens_R2.fastq molestus_R1.fastq molestus_R2.fastq )

SOURCEDIR='/data/joiglu/mosquito'
FILENAME=( Cpipiens_S1_L001_R1_001.fastq.gz
           Cpipiens_S1_L001_R2_001.fastq.gz
           Mol1-5_S1_L001_R1_001.fastq.gz
           Mol1-5_S1_L001_R2_001.fastq.gz )

# ========================================================================
# Below, there should be no absolute paths, which are computer-dependent.


DATADIR=`pwd | sed 's/results/data/'`
NEWNAME=( pip_R1 pip_R2 mol_R1 mol_R2 )
# NEWNAMEs must be in the same order as FILENAMEs.

if [ ! -d $DATADIR ]; then mkdir $DATADIR; fi
for i in 0 1 2 3; do
   # Note that the simplified dataset, used only for testing purposes, is not
   # compressed, while the real dataset is. We need to create symbolic links with
   # appropriate extensions:
   SUFIX='.fq'
   if echo ${FILENAME[$i]} | grep -q '.gz'; then SUFIX='.fq.gz'; fi
   if [ ! -e $DATADIR/${NEWNAME[$i]} ]; then
      ln -s $SOURCEDIR/${FILENAME[$i]} $DATADIR/${NEWNAME[$i]}$SUFIX
   fi
done

# Now, original fastq files are in $DATADIR, and they are named either pip_R1.fq, pip_R2.fq, etc.
# or pip_R1.fq.gz, pip_R2.fq.gz, etc., depending on whether they point to the simplified dataset
# or to the real one.


# I will create the barcode files for pipiens and molestus in this folder

if [ ! -e pipiens_barcode.txt ]; then
   echo -e 'GGTCGTAAATG\tPipFe1_R1.fastq\tPipFe1_R2.fastq'  > pipiens_barcode.txt
   echo -e 'CTAGTACCTG\tPipFe2_R1.fastq\tPipFe2_R2.fastq'  >> pipiens_barcode.txt
   echo -e 'AACTACGGG\tPipFe3_R1.fastq\tPipFe3_R2.fastq'    >> pipiens_barcode.txt
   echo -e 'TCGACGTT\tPipFe6_R1.fastq\tPipFe6_R2.fastq'    >> pipiens_barcode.txt
   echo -e 'GTCAGAGTATG\tPipMa4_R1.fastq\tPipMa4_R2.fastq' >> pipiens_barcode.txt
   echo -e 'ACTGAGACTG\tPipFe4_R1.fastq\tPipFe4_R2.fastq'  >> pipiens_barcode.txt
   echo -e 'CAGCTCTAG\tPipMa3_R1.fastq\tPipMa3_R2.fastq'   >> pipiens_barcode.txt
   echo -e 'TGATCTCG\tPipMa1_R1.fastq\tPipMa1_R2.fastq'    >> pipiens_barcode.txt
   echo -e 'AGCCATGAATG\tPipMa2_R1.fastq\tPipMa2_R2.fastq' >> pipiens_barcode.txt
   echo -e 'CCAATGCTTG\tPipMa5_R1.fastq\tPipMa5_R2.fastq'  >> pipiens_barcode.txt
   echo -e 'GATTGCAGG\tPipMa6_R1.fastq\tPipMa6_R2.fastq'   >> pipiens_barcode.txt
   echo -e 'TTGGCATC\tPipFe5_R1.fastq\tPipFe5_R2.fastq'    >> pipiens_barcode.txt
fi

if [ ! -e molestus_barcode.txt ]; then
   echo -e 'GGTCGTAAATG\tMol01_R1.fastq\tMol01_R2.fastq'  > molestus_barcode.txt
   echo -e 'GTCAGAGTATG\tMol02_R1.fastq\tMol02_R2.fastq' >> molestus_barcode.txt
   echo -e 'CTAGTACCTG\tMol03_R1.fastq\tMol03_R2.fastq'  >> molestus_barcode.txt
   echo -e 'AACTACGGG\tMol04_R1.fastq\tMol04_R2.fastq'   >> molestus_barcode.txt
   echo -e 'TCGACGTT\tMol05_R1.fastq\tMol05_R2.fastq'    >> molestus_barcode.txt
fi

if [ ! -e PipFe1_R1.fastq ]; then
   # I assume that sabre can take either compressed or uncompressed fastq files
   # as input
   sabre pe -m 1 -c \
            -f $DATADIR/pip_R1$SUFIX \
            -r $DATADIR/pip_R2$SUFIX  \
            -b pipiens_barcode.txt \
            -u pip_unknown_R1.fastq \
            -w pip_unknown_R2.fastq
fi

if [ ! -e Mol01_R1.fastq ]; then
   sabre pe -m 1 -c \
            -f $DATADIR/mol_R1$SUFIX \
            -r $DATADIR/mol_R2$SUFIX \
            -b molestus_barcode.txt \
            -u mol_unknown_R1.fastq \
            -w mol_unknown_R2.fastq
fi

LIST=(PipFe1 PipFe2 PipFe3 PipFe6 PipMa4 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5 Mol01 Mol02 Mol03 Mol04 Mol05)

# If pear is to run locally, it is better not to parallelize 17 processes, but 6
# at most. We just need to count how many processors there are:

PROC=`grep -P '^processor' /proc/cpuinfo | wc -l`

if [ ! -d merged ]; then mkdir merged; fi

# Now, I run on loop or the other, depending on whethere there are more or less
# than 8 processors available

if [ $PROC -gt 8 ]; then
   for i in `seq 0 16`; do
      if [ ! -e merged/${LIST[$i]}'_assembled.fastq' ]; then
         pear  -f ${LIST[$i]}'_R1.fastq' \
               -r ${LIST[$i]}'_R2.fastq' \
               -o merged/${LIST[$i]} \
               -v 10 \
               -q 15 \
               -j 1 \
               --memory 2G &
      fi
   done
   wait
else
   # With the loops below, no more than 6 processes will be running
   # at the same time
   for j in 0 6; do
      for i in `seq $j $(( j + 5 ))`; do
         if [ ! -e merged/${LIST[$i]}'_assembled.fastq' ]; then
            pear  -f ${LIST[$i]}'_R1.fastq' \
                  -r ${LIST[$i]}'_R2.fastq' \
                  -o merged/${LIST[$i]} \
                  -v 10 \
                  -q 15 \
                  -j 1 \
                  --memory 2G &
         fi
      done
      wait
   done
   for i in 12 13 14 15 16; do
      if [ ! -e merged/${LIST[$i]}'_assembled.fastq' ]; then
         pear  -f ${LIST[$i]}'_R1.fastq' \
               -r ${LIST[$i]}'_R2.fastq' \
               -o merged/${LIST[$i]} \
               -v 10 \
               -q 15 \
               -j 1 \
               --memory 2G &
      fi
   done
   wait
fi
