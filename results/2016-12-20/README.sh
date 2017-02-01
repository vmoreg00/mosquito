#!/bin/bash
#
#				2016-12-20
#				----------
#
# The goal of this analysis is to identify potential indel markers polymorphic
# pipiens and molestus. Previous analyses showed that individual samples of
# the pipiens form do not have much coverage per site. Thus, I resort to pool
# them. For that, I need to change their read groups and use freebayes again
# to identify variants between pipiens and molestus.

BAMDIR='../2016-07-12'
PICARDPATH=~/bin
SAMPLE=(Mol01 Mol02 Mol03 Mol04 Mol05 PipFe1 PipFe2 PipFe3 PipFe4 PipFe6 PipMa1 PipMa2 PipMa3 PipMa4 PipMa5 PipMa6)

if [ ! -e pip.vcf.gz ] || [ ! -e mol.vcf.gz ]; then
   if [ ! -e reference.fa ]; then
      wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz
      mv Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz reference.fa.gz
      gunzip reference.fa.gz
      samtools faidx reference.fa
   fi

   for i in "${SAMPLE[@]}"; do
      if [ ! -e $i.bam ]; then
         java -jar $PICARDPATH/picard.jar AddOrReplaceReadGroups \
         I=$BAMDIR/$i'_sorted.bam' \
         O=$i.bam \
         RGID=${i:0:3} \
         RGLB=${i:0:3} \
         RGPL=illumina \
         RGPU=MiSeq_IGB \
         RGSM=${i:0:3}
      fi
      if [ ! -e $i.bam.bai ]; then
         samtools index $i.bam
      fi
   done
   if [ ! -e pip.vcf.gz ]; then
      if [ ! -e pipbams.txt ]; then
         ls -1 Pip*.bam > pipbams.txt
      fi
      freebayes -f reference.fa \
                --bam-list pipbams.txt \
                --ploidy 11 \
                --pooled-discrete \
                --min-mapping-quality 15 \
                --min-base-quality 15 \
                --min-alternate-count 2 | \
      vcftools --vcf - --minQ 100 --recode --recode-INFO-all --stdout | \
      bgzip > pip.vcf.gz &
   fi
   if [ ! -e mol.vcf.gz ]; then
      if [ ! -e molbams.txt ]; then
         ls -1 Mol*.bam > molbams.txt
      fi
      freebayes -f reference.fa \
                --bam-list molbams.txt \
                --ploidy 5 \
                --pooled-discrete \
                --min-mapping-quality 15 \
                --min-base-quality 15 \
                --min-alternate-count 2 | \
      vcftools --vcf - --minQ 100 --recode --recode-INFO-all --stdout | \
      bgzip > mol.vcf.gz &
   fi
   wait
   tabix -p vcf pip.vcf.gz
   tabix -p vcf mol.vcf.gz
fi
