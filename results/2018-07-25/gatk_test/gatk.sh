#!/bin/bash

###############################################################################
#                                                                             #
# This script has been performed in order to check the GATK parameters with   #
# a subset of the original data                                               #
#                                                                             #
# Author: Moreno-Gonzalez, Victor                                             #
# Date:   2018-11-11                                                          #
###############################################################################

REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/data/asgharian/bam/bwa
GATK=/data/victor/src/GenomeAnalysisTK.jar  # version 3.8.1
sample=(1.molA1 2.pipA4 3.molM1 4.torM2 5.molS1 6.mixS2 7.mixS3 8.torM4) # Names
ploidy=(41 41 52 56 30 26 41 41); # Genome copies in each sample

# Subset
if [ ! -d bam ]; then
	mkdir bam;
fi;
for i in `seq 0 7`; do
	if [ ! -e bam/${sample[i]}.subset.bam ]; then
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

# Run GATK test
if [ ! -d vcf ]; then
        mkdir vcf;
fi;
# Detection of possible haplotypes
for i in `seq 0 7`; do
        if [ ! -e vcf/${sample[i]}.g.vcf ]; then
                java -jar $GATK -T HaplotypeCaller \
                        -R $REFGENOME/CulQui.fna \
                        -I bam/${sample[i]}.subset.bam \
                        --output_mode EMIT_ALL_CONFIDENT_SITES \
                        --emitRefConfidence GVCF \
                        --sample_ploidy ${ploidy[i]} \
                        --max_alternate_alleles 3 \
                        --max_num_PL_values 100000 \
                        --variant_index_type LINEAR \
                        --variant_index_parameter 128000 \
			--num_cpu_threads_per_data_thread 8 \
			--log_to_file vcf/${sample[i]}.log \
                        -o vcf/${sample[i]}.g.vcf &
        fi;
done;
wait

# Generation of the Genotypes and join the gvcf files
# (The Variant Recalibration cannot be done as SNPs are unknown in this species)
if [ ! -e vcf/culex_cohort.vcf ]; then
        java -Djava.io.tmpdir=tmp -Xmx60g -jar $GATK -T GenotypeGVCFs \
                -R $REFGENOME/CulQui.fna \
                --variant vcf/1.molA1.g.vcf \
                --variant vcf/2.pipA4.g.vcf \
                --variant vcf/3.molM1.g.vcf \
                --variant vcf/4.torM2.g.vcf \
                --variant vcf/5.molS1.g.vcf \
                --variant vcf/6.mixS2.g.vcf \
                --variant vcf/7.mixS3.g.vcf \
                --variant vcf/8.torM4.g.vcf \
		--num_threads 1 \
		--max_alternate_alleles 3 \
		--max_num_PL_values 1000000 \
		--annotateNDA \
		--log_to_file vcf/culex_cohort.log \
                -o vcf/culex_cohort.vcf;
fi;
