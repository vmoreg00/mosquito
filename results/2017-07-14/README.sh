#!/bin/bash
#
#				2017-07-14
#				----------
#
# On 2017-05-30 I downloaded the raw reads from some mosquitoes available
# and mapped them to the loci identified by ipyrad from our data. The next
# step  is to assemble the reads mapped into haplotypes. There are different
# options for this. An old post in a forum suggests GATK is better than samtools
# because of the better variation callers available in GATK. Here, I will use
# freebayes to call variants, hoping that freebayes will be able to keep good
# track of the two haplotypes potentially present in each individual. Plus,
# Michael G. Harvey has made available a script, freebayes_vcf2fa.py, that
# generates the fasta file from a vcf obtained with freebayes. By the way,
# I must provide this link to Harvey's repository in any subsequent publication:
#
# https://github.com/mgharvey/misc_python
#
# The problem is that the current version of freebayes, v1.1.0-44-gd784cf8,
# throws a "signal 11" error shortly after starting. An issue (#390) was raised
# in the github page, but not solved yet.

BAMDIR=../2017-05-30/bam
REFERENCE=../2017-05-30/consensus.fasta
ORIGINAL=( SRR2029627 SRR2029628 SRR2029629 SRR2029630
           SRR2029631 SRR2029632 SRR2029633 SRR2029634 )

NEWNAME=( PipA1 PipA4 PipM1 TorM2 PipS1  PipS2  PipS3 TorM4 )

for i in 0 1 2 3 4 5 6 7; do
   if [ ! -e ${NEWNAME[$i]}.fa ]; then
      if [ ! -e ${NEWNAME[$i]}.vcf ]; then
         freebayes -f $REFERENCE \
                   --report-monomorphic \
                   --max-complex-gap 50 \
                   $BAMDIR/${ORIGINAL[$i]}'_sorted.bam' > ${NEWNAME[$i]}.vcf
      fi
      python freebayes_vcf2fa.py ${NEWNAME[$i]}.vcf ${NEWNAME[$i]}.fa ${NEWNAME[$i]}.summary 2 1
   fi
done
