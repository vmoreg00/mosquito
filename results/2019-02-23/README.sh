###############################################################################
#                                                                             #
# ABBA/BABA test                                                              #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2019-02-23                                                            #
#                                                                             #
###############################################################################

VCF=/data/victor/mosquito/data/asgharian/vcf/FB_CRISP_filtered/culex_filtered.vcf

# Generation of the ABBA-BABA input files
if [ ! -e culex_alt_freqs.tsv ]; then
	python3 src/alt_allele_freqs.py $VCF \
                                        culex_alt_freqs.tsv \
                                        culex_chr_lengths.tsv;
fi;

# Genome-Wide ABBA-BABA test: 4 cases are tested
## Case 1: ((([mixS2+mixS3], molM1), pipA4), p0)
## Case 2: ((([mixS2+mixS3], molM1), molA1), p0)
## Case 3: (((molS1, molM1), pipA4), p0)
## Case 4: (((molS1, molM1), molA1), p0)
if [ ! -e abba_baba_1e+06.txt ]; then
	Rscript src/ABBA_BABA.R;
fi;

# Sliding-Windows ABBA-BABA test with non-overlaping 100kb windows
# Only cases 1 and 2 are tested.
if [ ! -e abba_baba_slidingWindows.png ]; then
	Rscript src/ABBA_BABA_SlidingWindows.R;
fi;

################################# Conclusions #################################
# In the genome-wide test, it is shown that when molestus is supposed to be
# the donor population (cases 4 and 2), the admixture proportion is higher
# than the cases where the donor is the pipiens form. This differences seem
# significant because the standar error of the statistic is very low.
#
# When looking to this differences in the sliding-windows analysis, the
# differences among cases 1 and 2 seems to be widespread distributed along
# the contigs, although in overall, the case 2 has higher f values than
# the case 1.
#
# The differences in the f-statistic between the cases, should be due to
# a preference in the mating by form.
