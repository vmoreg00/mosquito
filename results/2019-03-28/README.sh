###############################################################################
#                                                                             #
# Functional annotation and analysis of the ABBA/BABA test                    #
#                                                                             #
# Author: Moreno-Gonzalez, V.                                                 #
# Date: 2019-03-28                                                            #
#                                                                             #
###############################################################################

# WW test

# functional annotation
if [ ! -e diff_f-percentile99.tsv ]; then
	Rscript src/functional-annotation.R;
fi;
