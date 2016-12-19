#!/bin/bash
#
#				2016-07-25
#                               ----------
#
# Another way to analyze the alignment against the reference it is using a programme called Freebayes.
# It will allow us to do a complete test about the alignment, so we will be able to see SNPs and  variants genotypics.
# Calling population variants
DIR='../2016-07-12'
SAMPLE=(Mol01 Mol02 Mol03 Mol04 Mol05 PipFe1 PipFe2 PipFe3 PipFe4 PipFe5 PipFe6 PipMa1 PipMa2 PipMa3 PipMa4 PipMa5 PipMa6)

if [ ! -e fb.vcf.gz ]; then

   if [ ! -e bamlist.txt ]; then
      ls -1 $DIR/*.bam > bamlist.txt
   fi

   if [ ! -e populations.txt ]; then
      touch populations.txt
      for i in `seq 0 16`; do
         echo "${SAMPLE[$i]} ${SAMPLE[$i]:0:3}" >> populations.txt
      done
   fi

   if [ ! -e fb.vcf ]; then
      freebayes -f $DIR/reference.fa \
                -L bamlist.txt \
                -v fb.vcf \
                --min-alternate-count 1 \
                --populations populations.txt \
                --use-mapping-quality \
                --genotype-qualities \
                --min-mapping-quality 10 \
                --use-best-n-alleles 2
   fi

   bgzip -c fb.vcf > fb.vcf.gz
   tabix -p vcf fb.vcf.gz
   rm fb.vcf
   rm bamlist.txt
   # rm populations.txt
fi

#Run a test of Hardy-Weinberg equilibrium on each population separately.
if [ ! -e Pip.hwe ]; then
   vcftools --gzvcf fb.vcf.gz \
      --remove-indv Mol01 \
      --remove-indv Mol02 \
      --remove-indv Mol03 \
      --remove-indv Mol04 \
      --remove-indv Mol05 \
      --out Pip --hardy
fi

if [ ! -e Mol.hwe ]; then
   vcftools --gzvcf fb.vcf.gz \
      --indv Mol01 \
      --indv Mol02 \
      --indv Mol03 \
      --indv Mol04 \
      --indv Mol05 \
      --out Mol --hardy
fi

if [ ! -e summary_hwe.txt ]; then
   echo -e "Population\tGenotypes\tP_HWE\tFREQ" > summary_hwe.txt
   gawk '(NR > 1){F[$3 "\t" $6]++}END{for (f in F) print "pipiens \t" f "\t" F[f]}' Pip.hwe | \
   sort -nrk 4 >> summary_hwe.txt
   gawk '(NR > 1){F[$3 "\t" $6]++}END{for (f in F) print "molestus\t" f "\t" F[f]}' Mol.hwe | \
   sort -nrk 4 >> summary_hwe.txt
fi

if [ ! -e pipienspool.vcf.gz ]; then
   if [ ! -e pipiensbams.txt ]; then
      ls -1 $DIR/Pip*.bam > pipiensbams.txt
   fi
   if [ ! -e pipienspool.vcf ]; then
      freebayes -f $DIR/reference.fa \
                -L pipiensbams.txt \
                -p 24 \
                --pooled-discrete > pipienspool.vcf
   fi
   bgzip -c pipienspool.vcf > pipienspool.vcf.gz
   tabix -p vcf pipienspool.vcf.gz
   rm pipienspool.vcf
   rm pipiensbams.txt
fi
