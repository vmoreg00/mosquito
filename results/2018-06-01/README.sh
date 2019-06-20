#!/bin/bash

###############################################################################
#                                                                             #
# Once the Asgharian et al. 2015 data (Culex pipiens pipiens, C. pipiens      #
# molestus and C. torrentium from Sacramento (USA), Aleksin (Russia) and      #
# Moskow (Russia)) have been trimmed to ensure a convenient quality for the   #
# downstream analysis, I will map the sequences to the Culex quinquefasciatus #
# genome.                                                                     #
#                                                                             #
# To do this, I will download the reference genome from the NCBI, and then, I #
# will use Bowtie2 to map Asgharian et al. 2015 sequences.                    #
#                                                                             #
###############################################################################

DATA=/data/victor/mosquito/data/asgharian/trimmed
REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/data/asgharian/bam

EM=/data/victor/src/EM-SNP

#################  Downloading the C. quinquefasciatus genome #################
# C. quiquefasciatus genome has 579,042,118 nucleotides and 3,171 scaffolds.
# It was sequenced by Arensburger et al. 2010 and is supposed to have 3 pairs
# of chromosomes (Cq.1, Cq.2 and Cq.3), although it is still incompletely
# aligned.

if [ ! -d $REGENOME ]; then mkdir $REFGENOME; fi

if [ ! -e $REFGENOME/CulQui.fna ]; then
	# Download
	wget -o Cq.genome.log \
		ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/185/GCF_000209185.1_CulPip1.0/GCF_000209185.1_CulPip1.0_genomic.fna.gz \
		-P $REFGENOME;
	wget -o Cq.annot.log \
		ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/185/GCF_000209185.1_CulPip1.0/GCF_000209185.1_CulPip1.0_genomic.gff.gz \
		-P $REFGENOME;

	# Untar
	gunzip $REFGENOME/GCF_000209185.1_CulPip1.0_genomic.fna.gz;
	gunzip $REFGENOME/GCF_000209185.1_CulPip1.0_genomic.gff.gz;

	# Simplify the names
	ls -1 $REFGENOME | awk -v R=$REFGENOME '{print "mv "R"/"$0" "R"/"$0}' |
	sed 's/GCF_000209185.1_CulPip1.0_genomic/CulQui/2' | sh;
fi;


################ Mapping the C. pipiens and molestus sequences ################
# I have tried to map the sequences using both Bowtie2 and BWA to test which of
# them performs better.

if [ ! -d $BAM ]; then mkdir $BAM; fi;

# === bowtie2 test ===
# It is done in local mode to be less restrictive as I am trying to aligning
# different species. The similarity between them is spected high, but not very
# high.

## Index
if [ ! -e $REFGENOME/CulQui_index.1.bt2 ]; then
        bowtie2-build $REFGENOME/CulQui.fna \
                      $REFGENOME/CulQui_index \
                      --threads 30 \
                      1>index.log 2>index.err;
fi;
## Mapping
if [ ! -d $BAM/bt2 ]; then
	mkdir $BAM/bt2;
	mkdir bt2;
	for i in `seq 27 34`; do
		bowtie2 --local -N 1 --phred33 -p 30 \
			-x $REFGENOME/CulQui_index \
			-1 $DATA/SRR20296$i\_forward_paired.fq.gz \
			-2 $DATA/SRR20296$i\_reverse_paired.fq.gz \
			-S $BAM/bt2/SRR20296$i.sam > bt2/SRR20296$i.bowtie2.log
		samtools view -bS -o $BAM/bt2/SRR20296$i.bam $BAM/bt2/SRR20296$i.sam;
		rm $BAM/bt2/SRR20296$i.sam;
		samtools sort -o $BAM/bt2/SRR20296$i.sorted.bam $BAM/bt2/SRR20296$i.bam;
		samtools stats $BAM/bt2/SRR20296$i.sorted.bam > bt2/SRR20296$i.stats;
	done;
fi;

# === BWA test ===
## Index
if [ ! -d $REFGENOME/bwa ]; then
	mkdir $REFGENOME/bwa;
	bwa index -p $REFGENOME/bwa/CulQui_index $REFGENOME/CulQui.fna
fi;

## Mapping
if [ ! -d $BAM/bwa ]; then
	mkdir $BAM/bwa;
	mkdir bwa;
	for i in `seq 27 34`; do
		bwa mem $REFGENOME/bwa/CulQui_index \
			$DATA/SRR20296$i\_forward_paired.fq.gz \
			$DATA/SRR20296$i\_reverse_paired.fq.gz \
			-t 30 > $BAM/bwa/SRR20296$i.sam;
		samtools view -bS -o $BAM/bwa/SRR20296$i.bam $BAM/bwa/SRR20296$i.sam;
		rm $BAM/bwa/SRR20296$i.sam;
		samtools sort -o $BAM/bwa/SRR20296$i.sorted.bam $BAM/bwa/SRR20296$i.bam;
		samtools stats $BAM/bwa/SRR20296$i.sorted.bam > bwa/SRR20296$i.stats;
	done;
fi;

# After checking the *.stats files, I have decided to select the BWA mapping as
# the number of mapped reads and the coverage is greater in it.
