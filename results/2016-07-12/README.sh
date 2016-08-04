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
	samtools mpileup -gf reference.fa  *_sorted.bam > mosquito.bcf
	bcftools view -Acg mosquito.bcf > mosquito.vcf
#Run a test of Hardy-Weinberg equilibrium
if [ ! -e HW.hwe ]; then
   vcftools --vcf mosquito.vcf --out HW --hardy
fi

# There are 2498 sites where all samples are homozygous for an alternative
# SNP, and 1153 sites segregating among the samples. All the sites homozygous
# for the alternative allele carry the same information, namely:
#
# CHR     POS     OBS(HOM1/HET/HOM2)      E(HOM1/HET/HOM2)        ChiSq_HWE       P_HWE   P_HET_DEFICIT   P_HET_EXCESS
# supercont3.3138 668     0/0/17  0.00/0.00/17.00 -nan    1.000000e+00    1.000000e+00    1.000000e+00
#
# Because their position is not a concern right now, I erase the sites completely
# homozygous for the alternative allele. Not without counting the number of lines
# to erase. I intentionally make the execution of the following inconditional.

grep 0/0/17 HW.hwe | wc -l | tee NumSitesFixedAlt.txt
sed -i '/0\/0\/17/ d' HW.hwe
