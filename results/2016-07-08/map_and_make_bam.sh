#!/bin/bash
#
# syntax: map_and_make_bam.sh <fastq_dir> <sample_id> <clean_up>
#
# if a third argument is passed, intermediate files will be
# deleted, if present.

if [ ! -e $2'_r1r2.fastq' ]; then
   reformat.sh  in=$1/$2'_R1_trimmed.fastq' \
               in2=$1/$2'_R2_trimmed.fastq' \
               out=$2'_r1r2.fastq' \
               underscore
fi

if [ ! -e $2'.sam' ]; then
   bowtie2 --local \
           --very-sensitive \
           --rg-id $2 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$2 \
           -x culex \
           -U $2'_r1r2.fastq',$1/$2'_trimmed.fastq' \
           -S $2'.sam'
fi

samtools view -bS $2'.sam' 1> $2'.bam'

if [ $# -eq 3 ]; then
   if [ -e $2'_r1r2.fastq' ]; then rm $2'_r1r2.fastq'; fi
   if [ -e $2'.sam' ]; then rm $2'.sam'; fi
fi
