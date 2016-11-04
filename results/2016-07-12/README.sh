#!/bin/bash
#
#                               12-07-16
#                               ----------
#
#PATH
DIR=../2016-07-08
# For the visualization of SNPs and INDELS,
# we need to convert BAM files to VCF.
LISTA=(PipFe1 PipFe2 PipFe3 PipFe4 PipFe5 PipFe6 PipMa1 PipMa2 PipMa3 PipMa4 PipMa5 PipMa6 Mol01 Mol02 Mol03 Mol04 Mol05)

# Index the reference for SAMtools
# First, download the reference and then,
# if it does not exist, we index it by SAMtools faidx.
if [ ! -e reference.fa.fai ]; then
        wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz
	mv Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz reference.fa.gz
	gunzip reference.fa.gz
	samtools faidx reference.fa
fi
# Index and sort the reads.
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
	if [ ! -e  ${LISTA[$i]}'_sorted.bam' ]; then
		samtools sort  $DIR/${LISTA[$i]}.bam  ${LISTA[$i]}'_sorted'
	fi
	if [ ! -e ${LISTA[$i]}'_sorted.bam.bai' ]; then
		samtools index   ${LISTA[$i]}'_sorted'.bam
	fi
done
# Finally, we will be able to convert the BAM files obtained at the last step, to VCF.
# We will obtain a file with the information we called.
if [ ! -e mosquito.vcf.gz ]; then
	samtools mpileup -gDf reference.fa  *_sorted.bam > mosquito.bcf
	bcftools view -Acg mosquito.bcf > mosquito.vcf
	bgzip -c mosquito.vcf > mosquito.vcf.gz
	tabix -p vcf mosquito.vcf.gz
fi
#Run a test of Hardy-Weinberg equilibrium on each population separately.
if [ ! -e Pip.hwe ]; then
   vcftools --gzvcf mosquito.vcf.gz \
      --remove-indv Mol01 \
      --remove-indv Mol02 \
      --remove-indv Mol03 \
      --remove-indv Mol04 \
      --remove-indv Mol05 \
      --out Pip --hardy
fi

if [ ! -e Mol.hwe ]; then
   vcftools --gzvcf mosquito.vcf.gz \
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

