#!/bin/bash

###############################################################################
#                                                                             #
# When running freebayes, I indicate very few cores so the execution was too  #
# slow. After several weeks, I decided to stop the freebayes and re-run the   #
# remaining regions with more cores (60).                                     #
#                                                                             #
# This script performs the execution of freebayes with the remaining regions  #
# and concatenate the resulting VCF files.                                    #
#                                                                             #
# Author: Victor Moreno-Gonzalez                                              #
# Date: 2018-01-22                                                            #
#                                                                             #
###############################################################################

REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/data/asgharian/bam/bwa
VCF=/data/victor/mosquito/data/asgharian/vcf

FB=/data/victor/src/freebayes-parallel/freebayes-parallel

sample=(1.molA1 2.pipA4 3.molM1 4.torM2 5.molS1 6.mixS2 7.mixS3 8.torM4) # Names
ploidy=(41 41 52 56 30 26 41 41); # Genome copies in each sample

# Un-called regions
if [ ! -e regions_left.txt ]; then
	# Contruct the regions file
	python2 src/fasta_generate_regions.py $REFGENOME/CulQui.fna.fai 5000000 \
	       > tmp.txt;
	# Extract the cutting line
	LINES=`wc -l tmp.txt | cut -d " " -f 1`;
	REGION=`grep -n "NW_001887000" tmp.txt | cut -d ":" -f 1`;
	CUT=$((LINES - REGION + 1));
	# Extract the regions that are needed and remove tmp file
	tail -n $CUT tmp.txt > regions_left.txt;
	head -n $((REGION - 1)) tmp.txt > regions_done.txt;
	rm tmp.txt;
fi;

# Variant calling
if [ ! -e $VCF/FB/culex_cohort_2.vcf ]; then
        $FB regions_left.txt \
                60 -f $REFGENOME/CulQui.fna \
                --cnv-map $BAM/ploidies.txt \
                --pooled-discrete \
                --use-best-n-alleles 3 \
                --report-monomorphic \
                $BAM/culex_merged.bam > $VCF/FB/culex_cohort_2.vcf;
fi;

# Comparing and concatenating the VCFs
#
#   !!! DO NOT TO RUN !!!
#
# (Just for security, I prefer to run it manually)
#
# Compress and index original vcf files
if [ ! -e $VCF/FB/culex_cohort.vcf.gz ]; then
	bgzip -@ 50 $VCF/FB/culex_cohort.vcf;
	bcftools index --threads 10 $VCF/FB/culex_cohort.vcf.gz;
fi;

if [ ! -e $VCF/FB/culex_cohort_2.vcf.gz ]; then
	bgzip -@ 50 $VCF/FB/culex_cohort_2.vcf;
	bcftools index --threads 10 $VCF/FB/culex_cohort_2.vcf.gz;
fi;


# Modify the regions file for bcftools
if [ ! -e regions_done_bcftools.txt ]; then
	cat regions_done.txt \
	  | sed 's/:/\t/' \
	  | sed 's/-/\t/' \
	  > regions_done_bcftools.txt;
fi;

# Subset first 7000 contigs
if [ ! -e $VCF/FB/culex_cohort_filtered.vcf.gz ]; then
	bcftools view -o $VCF/FB/culex_cohort_filtered.vcf.gz -O z \
	              -R regions_done_bcftools.txt \
	              --threads 50 $VCF/FB/culex_cohort.vcf.gz;
fi;

# Concatenate
if [ ! -e $VCF/FB/culex_cohort_concat.vcf ]; then
	vcf-concat $VCF/FB/culex_cohort_filtered.vcf.gz \
		   $VCF/FB/culex_cohort_2.vcf.gz \
	           > $VCF/FB/culex_cohort_concat.vcf;
fi;

# Remove files (???)
rm regions*
