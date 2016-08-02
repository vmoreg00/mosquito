#!/bin/bash
#
# syntax: map_and_make_bam.sh <fastq_dir> <sample_id> <clean_up>
#
# if a third argument is passed, intermediate files will be
# deleted, if present.

if [ ! -e $2'.sam' ]; then
   bowtie2 --local \
           --very-sensitive \
           --rg-id $2 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$2 \
           -x culex \
           -U $1/$2'_setrimmed.fastq' \
           -S $2'.sam' &> $2.maplog
fi

samtools view -bS $2'.sam' 1> $2'.bam'

if [ $# -eq 3 ]; then
   if [ -e $2'.sam' ]; then rm $2'.sam'; fi
fi
