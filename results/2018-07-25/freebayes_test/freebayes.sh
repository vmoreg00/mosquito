#!/bin/bash

###############################################################################
#                                                                             #
# This script has been performe to check the freebayes parameters with a      #
# subset of the original data in order to get the desired output and a        #
# good performance                                                            #
#                                                                             #
# Author: Moreno-Gonzalez, Victor                                             #
# Date:   2018-11-13                                                          #
###############################################################################

REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/data/asgharian/bam/bwa
FB=/data/victor/src/freebayes-parallel/freebayes-parallel
REGIONS=/data/victor/src/freebayes-parallel/fasta_generate_regions.py
sample=(1.molA1 2.pipA4 3.molM1 4.torM2 5.molS1 6.mixS2 7.mixS3 8.torM4) # Names
ploidy=(41 41 52 56 30 26 41 41); # Genome copies in each sample

#========================== Subset the original data =========================#
if [ ! -d bam ]; then
        mkdir bam;
fi;
# Subset the BAM files (500,000 first reads)
# Creates also the cnv-map file of ploidies
for i in `seq 0 7`; do
        if [ ! -e bam/${sample[i]}.subset.bam ]; then
		# cnv-map
		echo ${sample[$i]} ${ploidy[$i]} >> ploidies.txt
		# Subset
                samtools view -H $BAM/${sample[i]}.sort.RG.markdup.bam \
                        > bam/${sample[i]}.subset.sam;
                samtools view $BAM/${sample[i]}.sort.RG.markdup.bam | \
                        head -n 1000000 >> bam/${sample[i]}.subset.sam;
                samtools view -bh bam/${sample[i]}.subset.sam \
                        > bam/${sample[i]}.subset.bam;
                samtools index -b bam/${sample[i]}.subset.bam;
                rm bam/${sample[i]}.subset.sam;
        fi &
done;
wait

# Merge BAMs
if [ ! -e bam/culex_merge.bam ]; then
	samtools merge bam/culex_merge.bam bam/*.bam;
	samtools index bam/culex_merge.bam bam/culex_merge.bam.bai -@ 10;
fi;

# To save time, I will select the needed regions in the ref-genome
if [ ! -e regions.txt ]; then
	# search tha last covered region in bam
	max_region=100000000;
	for i in `seq 0 7`; do
		last_region=`samtools view bam/${sample[$i]}.subset.bam | \
			tail -n 1`;
		last_region=`echo $last_region | \
			cut -d " " -f 3 | cut -d _ -f 2 | cut -d . -f 1`;
		if [ $last_region -lt $max_region ]; then
			max_region=$last_region;
		fi;
	done;
	# select only the covered regions
	line=`grep -n $max_region $REFGENOME/CulQui.fna.fai | cut -d : -f 1`;
	$REGIONS $REFGENOME/CulQui.fna.fai 30000000 > tmp.txt;
	head -n $line tmp.txt > regions.txt;
	rm tmp.txt
fi;

#============================= Variant Calling ===============================
if [ ! -d vcf ]; then
	mkdir vcf;
fi;

source activate py27;
if [ ! -e vcf/culex_merge.vcf ]; then
	$FB regions.txt 20 -f $REFGENOME/CulQui.fna \
		--cnv-map ploidies.txt \
		--pooled-discrete \
		--use-best-n-alleles 3 \
		--report-monomorphic \
		bam/culex_merge.bam > vcf/culex_merge.vcf;
fi;
source deactivate;
