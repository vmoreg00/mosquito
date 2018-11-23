#!/bin/bash

###############################################################################
#                                                                             #
# The Asgharian et al. 2015 data (Culex pipiens pipiens, C. pipiens molestus  #
# and C. torrentium from Sacramento (USA), Aleksin (Russia) and Moskow        #
# (Russia)) have been mapped against the reference genome of                  #
# C. quinquefasciatus using BWA as it performed better than Bowtie2.          #
#                                                                             #
# In this analysis the variant calling has been carried out with two callers: #
# GATK and freebayes. Finally, the outcomes of both programs has been merged  #
# assuming that the real detected SNPs are those detected by both algorithms. #
#                                                                             #
###############################################################################

REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/data/asgharian/bam/bwa
VCF=/data/victor/mosquito/data/asgharian/vcf

GATK=/usr/local/bin/GenomeAnalysisTK.jar
PICARD=/usr/local/bin/picard.jar

########################### Read groups assignation ###########################
# During the alignment I miss to add -R option in bwa mem, so now I have to add
# the read groups to be able to run the variant calling. I am using picard.jar
# for that purpose
for i in `seq 27 34`; do
	if [ ! -e $BAM/SRR20296$i\_RG.bam ]; then
		java -jar $PICARD AddOrReplaceReadGroups \
			I=$BAM/SRR20296$i.sorted.bam \
			O=tmp.bam \
			RGID=SRR20296$i \
			RGLB=SRR20296$i \
			RGPL=illumina \
			RGPU=unit1 \
			RGSM=SRR20296$i;
		samtools view -b -F 132 tmp.bam > $BAM/SRR20296$i\_RG.bam;
		rm tmp.bam;
	fi;
done;

############################## Variant Calling ################################
# For some reason, samtools generates empty .vcf files (only header). I will
# use GATK HaplotypeCaller and FreeBayes to detect SNPs as the coincident
# in both algorithms. FreeBayes has the advantage that can consider pooled
# sequencing.

# Dictionary creation
if [ ! -e $REFGENOME/CulQui.dict ]; then
	java -jar $PICARD CreateSequenceDictionary \
		R=$REFGENOME/CulQui.fna \
		O=$REFGENOME/CulQui.dict;
fi;

# Variant calling directory
if [ ! -d $VCF ]; then
	mkdir $VCF
fi;

## ... using GATK
if [ ! -d $VCF/GATK ]; then
	mkdir $VCF/GATK;
fi;
for i in `seq 27 34`; do
	if [ ! -e $VCF/GATK/SRR20296$i.vcf ]; then
		samtools index -b $BAM/SRR20296$i\_RG.bam
		java -jar $GATK -T HaplotypeCaller \
			-R $REFGENOME/CulQui.fna -I $BAM/SRR20296$i\_RG.bam \
			--emitRefConfidence GVCF \
			-variant_index_type LINEAR -variant_index_parameter 128000 \
			-o $VCF/GATK/SRR20296$i.vcf;
	fi;
done;

## ... using FreeBayes
#
# Following Ashgarian et al. (2015), the number of pooled individuals in each
# BAM file are:
# +-----------------------------------+----------------------------------+
# | ID          Sample   Individuals  | ID          Sample   Individuals |
# +-----------------------------------+----------------------------------+
# | SRR2029627  A1       224          | SRR2029631  S1       15          |
# | SRR2029628  A4       132          | SRR2029632  S2       13          |
# | SRR2029629  M1       26           | SRR2029633  S3       64          |
# | SRR2029630  M2       28           | SRR2029634  M4       195         |
# +-----------------------------------+----------------------------------+
#
# NOTE.1: Due to a memmory issue in some files (28, 33 and 34), the memmory
#         use exceded the RAM of the server and the process is automatically
#         killed without raising any error. I looked for a solution. According
#         with the FreeBayes GitHub site (https://github.com/ekg/freebayes/)
#         it seems that this is a very common problem when analysing pooled
#         samples and the --pooled-discrete option is activated. To solve this
#         issue, the author purpose the next:
#
#         "Users can now safely use --pooled-discrete provided they also use a
#         suitable setting for --use-best-n-alleles.  In practice, setting to 5
#         or lower should be sufficient to prevent memory blowup in most
#         situations.  For the time being, I suggest testing with progressively
#         lower settings, or simply setting it as low as you think reasonable.
#         (For pooled experiments focused on SNPs, this would be 2.)"
#
#         For that reasson, I use --use-best-n-alleles 3 to test if the memmory
#         issue is solved.
# NOTE.2: For samples SRR2029627, 28, 33 and 34 the last solution do not
#         seems to work. This is probably due to the high genome copies and
#         the low coverage in the alignment, what leads to an overflow in
#         memmory usage. Using freebayes-parallel do not solve the problem, so
#         I've decided to reduce the copy number of those samples to the
#         mean of the other ones (SRR2029629, 30, 31 and  32:
#         (26 + 28 + 15 + 13)*2 / 4 = 41 mean genome copies).

if [ ! -d $VCF/FB ]; then
	mkdir $VCF/FB;
fi;

# Culex is diploid, so in each sample it is expected the double of genome copies
ind=(41 41 52 56 30 26 41 41);      # Genome copies in each sample
j=0;                                # counter
for i in `seq 27 34`; do
	if [ ! -e $VCF/FB/SRR20296$i.vcf ]; then
		echo SRR20296$i
		freebayes -f $REFGENOME/CulQui.fna -p ${ind[$j]} \
			--pooled-discrete --use-best-n-alleles 3 \
			--report-monomorphic \
			$BAM/SRR20296$i\_RG.bam > $VCF/FB/SRR20296$i.vcf &
	fi;
	j=$((j+1));
done;
wait

############################## Quality control ################################
# With special atention in SNPs
# Output will be saved in 'stats' folder.
#./src/vcf_quality.sh

# RESULTS:
#
# Both algorithms has performed simirar variant calling in all the files.
# In general terms, the number of detected SNPs ranges from ~2 to ~4M and the
# ts/tv ratio, from ~1.3 to ~1.4. Note that the number of SNPs detected in the
# same sample by either algorithm are very similar, as well as the ts/tv.
# Similarly, the proportion of each subtitution type is very similar in
# the different samples and in both outcomes (see table bellow and pdf files in
# 'stats' folder).
#
# It is notorious the low number of SNPs detected in the sample SRR2029627.
# This sample will be discarded in the ongoing analysis.
#
# Focusing on the SNPs quality distribution with both algorithms, it is shown
# that, in general, this distribution is more smoothed in the freebayes'
# variant calling. Nevertheless, with this last algorithm there are more SNPs
# of QUAL == 0.
#
#                                                 +-------------------+-------------------+
#                                                 | GATK              | freebayes         |
# -----------+------------+------------+----------+-----------+-------+-----------+-------+
# Sample ID  | Species    | Location   | Habitat  | SNPs GATK | ts/tv | SNPs      | ts/tv |
# -----------+------------+------------+----------+-----------+-------+-----------+-------+
# SRR2029627 | molestus   | Aleksin    | Urban    | 130       | 1.16  | 276       | 1.31  |
# SRR2029628 | pipiens    | Aleksin    | Suburban | 2,959,732 | 1.38  | 3,135,539 | 1.35  |
# SRR2029629 | molestus   | Moskow     | Urban    | 4,602,923 | 1.34  | 4,050,541 | 1.34  |
# SRR2029630 | torrentium | Moscow     | Suburban | 3,472,638 | 1.37  | 3,199,071 | 1.35  |
# SRR2029631 | molestus   | Sacramento | Urban    | 2,506,767 | 1.35  | 2,617,778 | 1.34  |
# SRR2029632 | mixed      | Sacramento | Suburban | 4,343,097 | 1.39  | 3,931,428 | 1.42  |
# SRR2029633 | mixed      | Sacramento | Suburban | 3,890,560 | 1.33  | 3,904,363 | 1.32  |
# SRR2029634 | torrentium | Moskow     | Suburban | 2,000,416 | 1.40  | 2,000,739 | 1.36  |
# -----------+------------+------------+----------+-----------+-------+-----------+-------+


############################## SNPs filtering #################################
# The variants will be filter according to the following criteria:
#   - SNPs
#   - QUAL >= 30.0
#   - Coincidence between GATK and freebayes variant callings
# Sample SRR2029627 is removed from the analysis.
#
# The vcf filter that are going to be filtered are those obtained by freebayes
# as their information is more valuable and it quality is more homogeneus. In
# other words, the filter is going to be applied over the freebayes output.

#if [ ! -d $VCF/filt ]; then
#	mkdir $VCF/filt;
#fi;
#
#for i in `seq 28 34`; do
#	if [ ! -e $VCF/filt/SRR20296$i.vcf ]; then
#		vcftools --vcf $VCF/FB/SRR20296$i.vcf \
#			 --remove-indels --minQ 30.0 \
#			 --diff $VCF/GATK/SRR20296$i.vcf \
#			 --diff-site --out diff_SRR20296$i -c | \
#		awk -v OFS="\t" -v ORS="\n" \
#			'{ if ( NR != 1 ) { print $1, $2; } }' \
#			> $VCF/filt/SRR20296$i.pos;
#		src/filter_vcf_positions.py \
#			$VCF/FB/SRR20296$i.vcf \
#			$VCF/filt/SRR20296$i.pos \
#			$VCF/filt/SRR20296$i.vcf;
#	fi;
#done;


########################### Final Quality control #############################

#./src/vcf_quality.sh

# After filtering the sequences, the number of variants has been reduced to a
# 57-69% of the original ones. The Ts/Tv ration has raised in all the filtered
# samples and the substitution rates by base has been mantained. The SNPs
# quality distribution has homogenized and now there is any variant with
# QUAL < 30.0. The remaining variants are:
#
#                                      +-------------------+-------------------------+
#                                      | Before            | After                   |
# -----------+------------+------------+-----------+-------+-----------------+-------+
# Sample ID  | Species    | Location   | SNPs      | ts/tv | SNPs            | ts/tv |
# -----------+------------+------------+-----------+-------+-----------------+-------+
# SRR2029628 | pipiens    | Aleksin    | 3,135,539 | 1.35  | 1,793,440 (57%) | 1.43  |
# SRR2029629 | molestus   | Moskow     | 4,050,541 | 1.34  | 2,690,528 (66%) | 1.40  |
# SRR2029630 | torrentium | Moscow     | 3,199,071 | 1.35  | 2,143,240 (67%) | 1.45  |
# SRR2029631 | molestus   | Sacramento | 2,617,778 | 1.34  | 1,589,848 (61%) | 1.37  |
# SRR2029632 | mixed      | Sacramento | 3,931,428 | 1.42  | 2,285,119 (58%) | 1.44  |
# SRR2029633 | mixed      | Sacramento | 3,904,363 | 1.32  | 2,338,102 (60%) | 1.36  |
# SRR2029634 | torrentium | Moskow     | 2,000,739 | 1.36  | 1,373,283 (69%) | 1.48  |
# -----------+------------+------------+-----------+-------+-----------------+-------+
