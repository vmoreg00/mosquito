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

if [ ! -e variants.vcf ] && [ ! -e variants.vcf.gz ]; then
   $BCFTOOLS mpileup -f $REFERENCE --read-groups RG.txt --annotate FORMAT/AD -Ou $BAMDIR/*.bam | \
   $BCFTOOLS call --format-fields GQ,GP --variants-only --multiallelic-caller --prior 0.005 --output variants.vcf
fi

# From a fast look at a few loci in individual bam files, using samtools tview,
# I realize that some loci suffer from paralogy problems, but are not difficult to
# identify because they have a higher coverage and because of the presence of more
# than two haplotypes per individual. Unfortunately, mpileup does not make any effort
# to phase the variants. Up to now I've seen one program that can use reads to phase
# nearby variants, hapCUT2. But first, I need to filter the variants.
#
# A balanced allelic depth in heterozygous is a clear indication of a variant being
# true. If only homozygous samples are present, then a minimum depth should be required.
# The bcftools filter function offers some possibilities to use sample-specific fields.
# But I did not manage to sum the depths of reference alleles across samples, for example.
#
# Because each "contig" is a short locus, I'd rather filter first the loci than the
# SNPs themselves, since it is likely that all variants in a locus with paralogy are
# questionable. The histogram of average depth per locus (not shown) has a long tail
# of loci with very high coverage. Limiting the average depth at 500 removes 22 out
# 182 loci.
#
# This is an example of preamble line in the vcf file:

##contig=<ID=locus_203405,length=378>


if [ ! -e loci_dp_below_500.bed ]; then
   gawk '(/^##contig/){
      split($1,CHOPPEDLINE,/[=,>]/)
      LENGTH[CHOPPEDLINE[3]] = CHOPPEDLINE[5]
   }(/^[^#]/){
      split($8,INFO,/;/)
      for (i in INFO) {
         split(INFO[i],VAR,/=/)
         if (VAR[1] == "DP") {
            DEPTHSUM[$1] += VAR[2]
            NUMVAR[$1]++
         }
      }
   }END{
      for (LOCUS in DEPTHSUM) {
         if (DEPTHSUM[LOCUS] / NUMVAR[LOCUS] < 500) print LOCUS "\t" 0 "\t" LENGTH[LOCUS]
      }
   }' variants.vcf > loci_dp_below_500.bed
fi

if [ ! -e variants.vcf.gz ]; then
   bgzip variants.vcf
   tabix -p vcf variants.vcf.gz
fi

if [ ! -e filtered.vcf ]; then
   $BCFTOOLS filter --output filtered.vcf --regions-file loci_dp_below_500.bed  variants.vcf.gz
fi
