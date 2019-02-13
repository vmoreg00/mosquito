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
	# Bcftools requires a .gz file compresed with bgzip and indexed
	bgzip -@ 50 $VCF/culex_cohort_concat.vcf;
	bcftools index --threads 30 $VCF/culex_cohort_concat.vcf.gz;
	# Remove invariants
	bcftools view -o $VCF/culex_cohort_NoInvSites.vcf -O v \
                      --min-ac 1 --threads 50 \
		      $VCF/culex_cohort_concat.vcf.gz;
fi;

#=============================== Quality check ===============================
if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR;
fi;

# 1 - BCFTOOLS statistics
if [ ! -e $OUTDIR/freebayes_summary.pdf ]; then
        bcftools stats -F $REFGENOME/CulQui.fna \
                 $VCF/culex_cohort_NoInvSites.vcf > $OUTDIR/freebayes.stats;
        plot-vcfstats --prefix $OUTDIR/tmp/ \
                      --title freebayes \
                      --main-title "Freebayes Variant Calling" \
                      $OUTDIR/freebayes.stats;
        mv $OUTDIR/tmp/summary.pdf $OUTDIR/freebayes_summary.pdf
        rm -r $OUTDIR/tmp
fi;

# 2 - Quality distribution
if [ ! -e freebayes_quality_distribution.png ]; then
        python $PYSTATS $VCF/culex_cohort_NoInvSites.vcf \
	       $OUTDIR/freebayes_quality_distribution;
fi;

#================================ CONCLUSION ================================#
# There are 19,426,285 SNPs and 2,100,281 indels (21,526,566 SNVs in total).
# Among SNPs there are 6,329,629 multiallelic SNPs.
# The ts/tv ratio is 1.07 and the highest substitution rates are in the
# transitions C>T, G>A, A>G and T>C.
#
# Qualities reported by FreeBayes are high enough, althoug there are lots of
# SNVs with QUAL==0.
#
# The quality of the freebayes output is good and seems to report more
# multiallelic sites than CRISP. This should be a good reason to perform the
# ongoing analysis on the basis of freebayes' VCF, filtering it according
# to those sites that are in both VCF files.
