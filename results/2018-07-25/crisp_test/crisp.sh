# script to test the parameters in crisp and its possible parallelization
#
# author: Moreno-Gonzalez, V.
# date: 2018-12-24

# Input data: subsets used in ../gatk_test
REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/results/2018-07-25/gatk_test/bam # Small subsets

sample=(1.molA1 2.pipA4 3.molM1 4.torM2 5.molS1 6.mixS2 7.mixS3 8.torM4) # Names
ploidy=(41 41 52 56 30 26 41 41); # Genome copies in each sample

# file with pool size in each sample (2 x number of individuals, assuming diploid genomes)
for i in `seq 0 7`; do
	echo "$BAM/${sample[$i]}.subset.bam PS=${ploidy[$i]}" >> bam.list
done;

# variant calling
time crisp \
	--bams bam.list \
	--ref $REFGENOME/CulQui.fna \
	--VCF culex.vcf \
	--minc 1 \
	--mbq 5 \
	--mmq 10 \
	--perms 30000 \
	--filterreads 0 \
	--qvoffset 33 \
	--EM 1 \
	--verbose 2 > culex.log
