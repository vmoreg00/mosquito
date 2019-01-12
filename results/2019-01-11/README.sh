#!/bin/bash

##############################################################################
#                                                                            #
# Variant calling performed by freebayes is analysed. As I indicated         #
# (perhaps erronously) that the program must return the invariable loci, the #
# final VCF is very heavy.                                                   #
#                                                                            #
# Author: Moreno-Gonzalez, V.                                                #
# Date: 2019-01-11                                                           #
#                                                                            #
##############################################################################

REFGENOME=/data/victor/mosquito/data/refgenome
VCF=/data/victor/mosquito/data/asgharian/vcf/FB
OUTDIR=/data/victor/mosquito/results/2019-01-11/stats
PYSTATS=/data/victor/mosquito/results/2019-01-07/src/snp_qual_dist.py

#=============================== Preprocessing ===============================
# As the VCF is very heavy, first of all I have to remove the invariable sites
# in order to handle this file more quickly.

# 1 - Remove non-variant sites
if [ ! -e $VCF/culex_cohort_NoInvSites.vcf ]; then
	bcftools view $VCF/culex_cohort.vcf \
                 -o $VCF/culex_cohort_NoInvSites.vcf -O v \
                 --min-ac 1 --threads 1;
fi;

# 2 - Compress original VCF file
if [ -e $VCF/culex_cohort.vcf ]; then
	bcftools view $VCF/culex_cohort.vcf -o $VCF/culex_cohort.bcf -O b;
	rm $VCF/culex_cohort.vcf;
fi;

#=============================== Quality check ===============================
if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR;
fi;

# 1 - BCFTOOLS statistics
if [ ! -e $OUTDIR/freebayes_summary.pdf ]; then
        bcftools stats -F $REFGENOME/CulQui.fna \
                 $VCF > $OUTDIR/freebayes.stats;
        plot-vcfstats --prefix $OUTDIR/tmp/ \
                      --title freebayes \
                      --main-title "Freebayes Variant Calling" \
                      $OUTDIR/freebayes.stats;
        mv $OUTDIR/tmp/summary.pdf $OUTDIR/freebayes_summary.pdf
        rm -r $OUTDIR/tmp
fi;

# 2 - Quality distribution
if [ ! -e freebayes_quality_distribution.png ]; then
        python $PYSTATS $VCF $OUTDIR/freebayes_quality_distribution;
fi;

