#!/bin/bash
#
#				2016-07-25
# ----------
# Another way to analyze the alignment against the reference it is using a programme called Freebayes.
# It will allow us to do a complete test about the alignment, so we will be able to see SNPs and  variants genotypics.
# Calling population variants
DIR='../2016-07-12'
if [ ! -e archivos.txt ]; then
	ls -1 ../2016-07-12/*.bam > archivos.txt
fi
if [ ! -e test.vcf ]; then
	freebayes  -f $DIR/reference.fa -L archivos.txt -v test.vcf  --use-mapping-quality --genotype-qualities --min-mapping-quality 20 -n 2
fi

# Run a test of Hardy-Weinberg equilibrium
vcftools --vcf test.vcf --out HW --hardy

