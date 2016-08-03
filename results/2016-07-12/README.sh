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
	if [ ! -e  ${LISTA[$i]}'_sorted' ]; then
		samtools sort  $DIR/${LISTA[$i]}.bam  ${LISTA[$i]}'_sorted'
	fi
	if [ ! -e ${LISTA[$i]}'_sorted'.bam.bai ]; then
		samtools index   ${LISTA[$i]}'_sorted'.bam
	fi
done
# Finally, we will be able to convert the BAM files obtained at the last step, to VCF.
# We will obtain a file with the information we called.
if [ ! -e mosquito.vcf ]; then
	samtools mpileup -gf reference.fa  *_sorted.bam > mosquito.mpileup
	bcftools view -Acg mosquito.mpileup > mosquito.vcf
fi
#Run a test of Hardy-Weinberg equilibrium 
vcftools --vcf mosquito.vcf --out HW --hardy

