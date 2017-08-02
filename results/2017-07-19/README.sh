#!/bin/bash
#
#				2017-07-19
#				----------
#
# On 2017-05-30 I downloaded the raw reads from some mosquitoes available
# and mapped them to the loci identified by ipyrad from our data. The next
# step  is to assemble the reads mapped into haplotypes. There are different
# options for this. An old post in a forum suggests GATK is better than samtools
# because of the better variation callers available in GATK. Here, I will use
# bcftools.

BAMDIR=../2017-05-30/bam
REFERENCE=../2017-05-30/consensus.fasta
BCFTOOLS=/usr/local/bin/bcftools
ORIGINAL=( SRR2029627 SRR2029628 SRR2029629 SRR2029630
           SRR2029631 SRR2029632 SRR2029633 SRR2029634 )

NEWNAME=( PipA1 PipA4 PipM1 TorM2 PipS1  PipS2  PipS3 TorM4 )

# It would make sense to call variants within each species, Culex pipiens and Culex
# torrentium, if the algorithm used a population model. Bcftools call can take
# advantage of known allele frequencies, but we can hardly estimate them on the fly
# with only a handful of individuals, anyways. So, I will call variants in all
# individuals together.

if [ ! -e RG.txt ]; then
   for i in 0 1 2 3 4 5 6 7; do
      echo "${ORIGINAL[$i]} ${NEWNAME[$i]}" >> RG.txt
   done
fi

if [ ! -e variants.vcf ]; then
   $BCFTOOLS mpileup -f $REFERENCE --read-groups RG.txt --annotate FORMAT/AD -Ou $BAMDIR/*.bam | \
   $BCFTOOLS call --format-fields GQ,GP --variants-only --multiallelic-caller --prior 0.005 --output variants.vcf
fi
