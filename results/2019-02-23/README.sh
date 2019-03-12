###############################################################################
#                                                                             #
# ABBA/BABA test                                                              #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2019-02-23                                                            #
#                                                                             #
###############################################################################

VCF=/data/victor/mosquito/data/asgharian/vcf/FB_CRISP_filtered/culex_filtered.vcf

if [ ! -e culex_alt_freqs.tsv ]; then
	python3 src/alt_allele_freqs.py $VCF \
                                        culex_alt_freqs.tsv \
                                        culex_chr_lengths.tsv;
fi;

if [ ! -e abba_baba_1e+06.txt ]; then
	Rscript src/ABBA_BABA.R;
fi;
