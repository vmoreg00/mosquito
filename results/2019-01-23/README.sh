#!/bin/bash

###############################################################################
#                                                                             #
# After the variant calling with FreeBayes and CRISP, it is necessary to      #
# filter the VCFs. This will be achieved by considering the true SNPs as      #
# those that have been called by both programs and that have a minimum        #
# quality of 15 in FreeBayes' VCF.                                            #
#                                                                             #
# After the filtering, the TreeMix's input file is made from the filtered VCF #
# file. Cx. torrentium populations will be treated as outgroup in the ML      #
# graph.                                                                      #
#                                                                             #
# Author: Victor Moreno-Gonzalez                                              #
# Date: 2018-01-23                                                            #
#                                                                             #
###############################################################################

VCF=/data/victor/mosquito/data/asgharian/vcf

FILTER=/data/victor/mosquito/results/2018-07-25/src/filter_vcf_positions.py

#======================= Filtering VCFs by difference =======================#
# Filtering the subset with CRISP's VCF file and by quality (minQ == 15)
if [ ! -d $VCF/FB_CRISP_filtered ]; then
	mkdir $VCF/FB_CRISP_filtered;
fi;

if [ ! -e $VCF/FB_CRISP_filtered/culex_filtered.vcf ]; then
	# CRISP have not find SNPs in all chromosomes, so there are some of
        # them that are absent. vcftools will fail if chromosomes are not the
        # same in both files. Thus, I have to search which chromosomes are
        # in CRISP VCF file and indicate them in vcftools execution with
        # --chr option.

	# Chromosomes in CRIPS
	grep -a -v "#" $VCF/crisp/culex_cohort.vcf | grep -a "^NW" | cut -f 1 | uniq \
	  > chr_in_crisp.txt;

	# Differences
	# --remove-indels options is not idicated because it is problematic.
	vcftools --vcf $VCF/FB/culex_cohort_NoInvSites.vcf \
	         `awk -v ORS=" " '{ print " --chr "$0 }' chr_in_crisp.txt` \
	         --minQ 15 \
	         --diff $VCF/crisp/culex_cohort.vcf \
	         --diff-site --out $VCF/FB_CRISP_filter/culex_diff -c \
	  | awk -v OFS="\t" -v ORS="\n" \
	        '{ if ( NR != 1 && $2 != "." && $3 != "." ) { print $1, $2 } }' \
          > culex_FB-crisp_common-possitions.pos;
	python $FILTER $VCF/FB/culex_cohort_NoInvSites.vcf \
	               culex_FB-crisp_common-possitions.pos \
	               $VCF/FB_CRISP_filtered/culex_filtered.vcf;
fi;

#============================= TreeMix input file ============================#
if [ ! -e treemix_input.tsv.gz ]; then
	grep -v "^#" $VCF/FB_CRISP_filtered/culex_filtered.vcf \
	  | cut -f 10-17 > culex_genotypes.tsv;

	src/treemix_input_creator.py \
		culex_genotypes.tsv \
		treemix_input.tsv;
	gzip treemix_input.tsv;
        gzip culex_genotypes.tsv;
fi;

#============================= TreeMix execution =============================#
if [ ! -d culex_treemix_out ]; then
	mkdir culex_treemix_out;
fi;
# Maximum Likelihood tree
if [ ! -e culex_treemix_out/culex_stem.treeout.gz ]; then
	treemix -i treemix_input.tsv.gz \
	        -o culex_treemix_out/culex_stem;
fi;

# Visualization
if [ ! -e culex_treemix_out/culex_tree.pdf ]; then
        # poporder file
        pops=(1.molA1 2.pipA4 3.molM1 5.molS1 6.mixS2 7.mixS3 4.torM2 8.torM4)
        for i in `seq 0 7`; do
                echo ${pops[i]} >> poporder;
        done;
        # plot phylogeny and residues
	Rscript src/treemix_plot.R;
        rm poporder;
fi;

# Bootstrap support
if [ ! -d culex_treemix_bootstrap ]; then
        mkdir culex_treemix_bootstrap;
fi;
if [ ! -e culex_treemix_bootstrap/culex_bootstrap_support.pdf ]; then
	# Bootstrap
	for i in {1..100}; do
	        treemix -i treemix_input.tsv.gz \
			-bootstrap -k 500 \
	                -o culex_treemix_bootstrap/culex_stem;
		mv culex_treemix_bootstrap/culex_stem.treeout.gz \
		   culex_treemix_bootstrap/bootstrap_$i.gz;
	done;
	rm culex_treemix_bootstrap/culex_stem.*;
	# Plot
        pops=(1.molA1 2.pipA4 3.molM1 5.molS1 6.mixS2 7.mixS3 4.torM2 8.torM4)
        for i in `seq 0 7`; do
                cat ${pops[i]} >> poporder;
        done;
        Rscript src/treemix_plot_bootstrap.R;
        rm poporder;
fi;

#================================= Conclusion ================================#
# After the excution of treemix, it is stated that the considered populations
# has tree-like historical relations and no migration events has occured.
# The tree topology shows that the Culex pipiens forms are more related by
# their location instead of their ecotype. This can indicate that the ecotype
# has no taxonomical support, leading to the conclussion that these 'forms' are
# just different phenotypes.
