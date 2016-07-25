#!/bin/bash
#
#				2016-07-25
# ----------
# Another way to analyze the alignment against the reference it is using a programme called Freebayes.
# It will allow us to do a complete test about the alignment, so we will be able to see SNPs and  variants genotypics.
# Calling population variants: Add read groups to BAM input
#Install:git clone --recursive https://github.com/ekg/bamaddrg.git
#cd bamaddrg
#make
DIR='../2016-07-12'
SUFIX='_sorted.bam.bai'

# First of all, we need to add read groups by bamaddrg
# after that, we execute freebayes where we need 2 types of input: reference and alignment bam file.
if [ ! -e test.vcf ]; then
	bamaddrg -b $DIR/PipFe1$SUFIX -s pip -r g1.1 \
               -b $DIR/PipFe2$SUFIX -s pip-r g1.2 \
               -b $DIR/PipFe3$SUFIX  -s pip -r g1.3 \
               -b $DIR/PipFe4$SUFIX  -s pip -r g1.4 \
	       -b $DIR/PipFe5$SUFIX  -s pip -r g1.5 \ 
	       -b $DIR/PipFe6$SUFIX -s pip -r g1.6 \
	       -b $DIR/PipMa1$SUFIX -s pip -r g1.7 \
	       -b $DIR/PipMa2$SUFIX -s pip -r g1.8 \
	       -b $DIR/PipMa3$SUFIX -s pip -r g1.9 \
	       -b $DIR/PipMa4$SUFIX -s pip -r g1.10 \
	       -b $DIR/PipMa5$SUFIX -s pip -r g1.11 \
	       -b $DIR/PipMa6$SUFIX -s pip -r g1.12 \
	       -b $DIR/Mol01$SUFIX -s mol\ -r g2.1 \
	       -b $DIR/Mol02$SUFIX -s mol\ -r g2.2 \
               -b $DIR/Mol03$SUFIX -s mol\ -r g2.3 \
	       -b $DIR/Mol04$SUFIX -s mol\ -r g2.4 \
               -b $DIR/Mol05$SUFIX -s mol\ -r g2.5 \
	 | freebayes -f $DIR/reference.fa --bam  $DIR/PipFe2$SUFIX $DIR/PipFe3$SUFIX $DIR/PipFe4$SUFIX  $DIR/PipFe5$SUFIX  $DIR/PipFe6$SUFIX $DIR/PipMa1$SUFIX $DIR/PipFe1$SUFIX \
	$DIR/PipMa2$SUFIX $DIR/PipMa3$SUFIX $DIR/PipMa4$SUFIX  $DIR/PipMa5$SUFIX  $DIR/PipMa6$SUFIX  $DIR/Mol01$SUFIX $DIR/Mol02$SUFIX  $DIR/Mol03$SUFIX $DIR/Mol04$SUFIX \
	$DIR/Mol05$SUFIX -v test.vcf --min-coverage 5 --use-mapping-quality --genotype-qualities
fi
