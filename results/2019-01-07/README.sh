#!/bin/bash

##############################################################################
#                                                                            #
# After the variant calling performed with CRISP I am studying the quality   #
# distribution of each variant.                                              #
#                                                                            #
# According to the program's log, CRISP has detected 15,666,187 SNVs in the  #
# 8 samples. It has token 8d 14H 49min 55.5s.                                #
#                                                                            #
# Author: Moreno-Gonzalez, V.                                                #
# Date: 2019-01-07                                                           #
#                                                                            #
##############################################################################

REFGENOME=/data/victor/mosquito/data/refgenome
VCF=/data/victor/mosquito/data/asgharian/vcf/crisp/culex_cohort.vcf
OUTDIR=/data/victor/mosquito/results/2019-01-07/stats

if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR;
fi;

# BCFTOOLS statistics
if [ ! -e $OUTDIR/crisp_summary.pdf ]; then
        bcftools stats -F $REFGENOME/CulQui.fna \
                 $VCF > $OUTDIR/crisp.stats;
        plot-vcfstats --prefix $OUTDIR/tmp/ \
                      --title crisp --main-title "CRISP Variant Calling" \
                      $OUTDIR/crisp.stats;
	mv $OUTDIR/tmp/summary.pdf $OUTDIR/crisp_summary.pdf
	rm -r $OUTDIR/tmp
fi;

# Generates the Quality distribution
if [ ! -e crisp_quality_distribution.png ]; then
	python src/snp_qual_dist.py $VCF ./stats/crisp_quality_distribution;
fi;

#================================ CONCLUSION ================================#
# There are 13,885,349 of SNPs and 2,100,281 indels (15,985,630 SNVs in
# total). The reassons why this number differ from the one given by the
# CRISP's log are unknown.
# Among the SNPs there are only 527,499 multiallelic SNPs.
# The ts/tv ratio is 1.23 and the highest substitution rates are in the
# transitions C>T, G>A, A>C and T>C.
#
# The qualities reported by CRISP are very high. The minimum reported
# quality is 20 while the maximum is 4,045. The mean quality is arround 100.
#
# For detailed information and figures, check 'stats' folder.
