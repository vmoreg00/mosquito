#!/bin/bash

###############################################################################
#                                                                             #
# In this analysis, I carry out a Fst analysis as the one preformed by        #
# Asgharian et al. (2015). The aims of this analysis are: (1) check the       #
# validity of the historical relationships obtained by TreeMix, and           #
# (2) compare my results with those obtained by Asgharian et al. (2015) in    #
# their original analysis.                                                    #
#                                                                             #
# I had have to develop a python script to transform the VCF file to the      #
# input format required by PoPoolation2, as the scripts provided by the       #
# developers do not accept VCF files with multiple populations.               #
#                                                                             #
# When computing the average Fst across all the genome, I have discarded      #
# all Fst values of those windows with less than 10 SNPs                      #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2019-06-09                                                            #
#                                                                             #
###############################################################################

VCF=/data/victor/mosquito/data/asgharian/vcf/FB_CRISP_filtered/culex_filtered.vcf
fst=/data/victor/src/popoolation2_1201/fst-sliding.pl
sync=src/multivcf2sync.py
avgFst=src/avg-fst.py
plot=src/fst_plot.R

# Transform form VCF to SYNCHRONIZE format
if [ ! -e culex.sync ]; then
	python $sync $VCF culex.sync;
fi;

# Compute the Fst statistic in 10Kb non-overlapping sliding-windows
if [ ! -e culex.fst ]; then
	./popoolation2_1201/fst-sliding.pl \
		--input culex.sync \
		--output culex.fst \
		--min-coverage 2 \
		--max-coverage 500 \
		--window-size 10000 \
	        --step-size 10000 \
	        --pool-size 41:41:52:56:30:26:41:41;
fi;

# Compute mean and standard deviation across the whole genome
if [ ! -e culex_fst_mean.csv ]; then
	python $avgFst culex.fst culex_fst
fi;

# Plot the results
if [ ! -e culex_fst_phylo.png ]; then
	Rscript $plot
fi;

#================================ CONCLUSION =================================#
# In a global overview, the results are quite similar to those obtained with
# TreeMix, although a small differences has been found in the "pipiens" clade.
# In the TreeMix graph, the sample molM1 is closer to the Sacramento samples,
# while in this new analysis, the samples from Alexin are the closest to
# the Sacramento samples.
#
# The mean Fst values are very similar to those obtained by Asgharian et al.
# (2015). Although my results show higher Fst values, in proportion, all
# Fst comparissons has the same proportion than in asgharian, so the
# relations in the dendrogram are the same.
#=============================================================================#
