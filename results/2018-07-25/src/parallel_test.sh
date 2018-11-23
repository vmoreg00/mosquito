##############################################################################
#                                                                            #
# parallel_test.sh                                                           #
#                                                                            #
# Script to test the performance of parallel-freebayes with a small bam      #
#                                                                            #
# File to test: SRR2029628.bam                                               #
#                                                                            #
# Author: Victor Moreno-Gonzalez                                             #
# Date: 09-10-2018                                                           #
#                                                                            #
##############################################################################

# When running FreeBayes, some processes were killed when running the SNP
# detection on files SRR2029627, 28 and 34 respectively. The reason may be a
# memmory issue that cannot be fixed with the parameter --use-best-n-alleles.
#
# The reference genome of Culex quinquefasciatus has a length of 586632570 nt.
# Splitting it into 5,000,000 nt length fragments and running 5 processes in
# parallel leads to 24 executions.
#
# I will run the test with SRR2029628.vcf as it was one of the files that
# fail to run freebayes.

VCF=/data/victor/mosquito/results/2018-07-25/prueba/SRR2029628.vcf
BAM=/data/victor/mosquito/data/asgharian/bam/bwa/SRR2029628_RG.bam
REFGENOME=/data/victor/mosquito/data/refgenome/CulQui.fna

if [ ! -d src ]; then
	mkdir src;
fi;

# Download the scripts of freebayes parallel
## NOTE: freebayes-parallel have been modifyied to fit the path to the other
##       scripts.
if [ ! -e src/freebayes-parallel ]; then
	wget https://raw.githubusercontent.com/ekg/freebayes/master/scripts/freebayes-parallel \
		-O src/freebayes-parallel;
	wget https://raw.githubusercontent.com/ekg/freebayes/master/scripts/fasta_generate_regions.py \
		-O src/fasta_generate_regions.py;
	wget https://raw.githubusercontent.com/vcflib/vcflib/5e3ce04f758c6df16bc4d242b18a24d725d2e6e5/scripts/vcffirstheader
		-O src/vcffirstheader
	chmod +x src/*
fi;

# Run fasta_generate_regions.py to see its output
if [ ! -e prueba/test.regions.fai ]; then
	src/fasta_generate_regions.py $REFGENOME.fai 5000000 > prueba/test.regions.fai;
fi;

# Run freeyes parallel (5 threads; 5000000 nt regions)
if [ ! -e $VCF ]; then
	time \
	src/freebayes-parallel <(src/fasta_generate_regions.py $REFGENOME.fai 5000000) \
			5 -f $REFGENOME -p 128 \
			--pooled-discrete --use-best-n-alleles 3 \
			$BAM > $VCF;
fi;

################################## Conclusions ################################
#
# After more than 8 days (11922m27.003s), freebayes parallel only had covered
# about 2000 regions of the 3069. As each day was slower, I've decided to stop
# it due to the expended time. If each file has to run alogn more than 20 days,
# as I'd had to run it with 4 files (see above), it will suppose nearly 80 days.
#
# Taking in account that this waiting time is not reasonable, I have decided to
# run freebayes again, but indicating that there are fewer genomic copies and
# selecting only the 3 best alleles.
#
