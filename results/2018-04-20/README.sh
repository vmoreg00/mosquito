###############################################################################
# The following are the accession numbers of the raw read from the study by
# Asgharian et al. 2015, Evolutionary genomics of Culex pipiens: global and
# local adaptations associated with climate, life-history traits and
# anthropogenic factors; Proc. R. Soc. B. 282(1810).
# Here there is some information about the reads:
#
# -------------------------------------------------------------------------------
# Run           Sample	Organism        	Sampling_site           MegaBases
# -------------------------------------------------------------------------------
# SRR2029627	A1	C. pipiens f. molestus	Aleksin Urban		2516
# SRR2029628	A4	C. pipiens f. pipiens	Aleksin Suburban	4901
# SRR2029629	M1	C. pipiens f. molestus	Moscow Urban		6618
# SRR2029630	M2	C. torrentium		Moscow Suburban		6535
# SRR2029631	S1	C. pipiens f. molesuts	Sacramento Urban	3370
# SRR2029632	S2	C. pipiens both forms	Sacramento Suburban	6037
# SRR2029633	S3	C. pipiens both forms	Sacramento Suburban	5520
# SRR2029634	M4	C. torrentium		Moscow Suburban		3306
# -------------------------------------------------------------------------------
#
# There is hardly any information on how the authors processed the reads.
# It may be necessary to trim and filter the reads before mapping them.
#
# In my opinion, as the aim of the job is to determine wheter there has been
# introgression or not between pipiens and molestus forms, it is not necessary
# to include the sequences of C. torrentium as it do not belong to the
# "Pipiens" complex (it is formed by C. pipiens (f. pipiens and f. molestus),
# C.quinquefasciatus, and, possibly, C. australicus. (see Harbach, 2012;
# doi.org/10.2987/8756-971X-28.4.10).
#
##############################################################################



###################### Downloading and checking quality ######################

# Downloading and looking the quality of the raw reads of Asgharian et al.
# Be patient because it will take arround 12 hours!
# Needed: sra-tools and fastqc
# Probably, the data of SRR2029632 (S2) and SRR2029633 (S3) will not be used
# because they belong to pooled samples of molestus and pipiens forms.

DATA=/data/victor/mosquito/data/asgharian/raw
SRATOOL=/data/victor/src/sratoolkit.2.9.0/bin

if [ ! -d $DATA ]; then mkdir $DATA; fi
if [ ! -d $DATA/fastqc ]; then mkdir $DATA/fastqc; fi

for i in SRR2029627 SRR2029628 SRR2029629 SRR2029630 \
         SRR2029631 SRR2029632 SRR2029633 SRR2029634; do
   if [ ! -e $DATA/$i\_1.fastq ]; then
      $SRATOOL/fastq-dump -F --split-files --outdir $DATA $i
   fi
   if [ ! -e $DATA/fastqc/$i\_1_fastqc.zip ]; then
      fastqc $DATA/$i\_1.fastq --outdir ./fastqc &
   fi
   if [ ! -e $DATA/fastqc/$i\_2_fastqc.zip ]; then
      fastqc $DATA/$i\_2.fastq --outdir ./fastqc &
   fi
done

# The files contains 392,148,848 reads (although Asgharian et al. say that
# there are 407 millions...). They were "sequenced as paired-ended 101bp
# reads on an Illumina HiSeq".

# After checking for the quality, it has been noted that, in general, there is
# a quality fall in the last ~ 10 bases of the reads. In adition, there are
# presence of overrepresented sequences, some of them belongs to Illumina
# adapters. In SRR2029628 (A4), there have also been found 9 overrepresented
# sequences with no hit and are characterized by short tandem repeats. In fact,
# the 9 sequences correspond with the repetition of AAGAATTT, AAAAATAC and
# CTTAAATT. In total, they correspond with approximately the 5% of all the
# reads in SRR2029628.
#
# ==> Found Adapters:
# TrueSeq Adapter (begening) 5'           GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
# Illumina Multiplexin PCR Primer 2.01    GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
# Illumina Single End PCR Primer 1        AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#
# ==> Other overrepresented sequences:
# AAGAATTT
# AAAAATAC
# CTTAAATT
#
# To remove all the overrepresented adapters, I will use trimmomatic, which
# trims the Illumina adapters, trims the extreme low quality bases and removes
# the low quality reads. For now, I will keep the unidentifyied overrepresented
# sequences.


################################## Trimming ###################################

# Trimmomatic pipeline: With trimmomatic, I will trim the adapters and the
# poor quality bases and reads. After, that, I will run fastqc again to check
# if the quality has improved.

ADAPTERS=/data/victor/src/Trimmomatic-0.38/adapters
FASTQ=/data/victor/mosquito/data/asgharian/raw
TRIMMED=/data/victor/mosquito/data/asgharian/trimmed
TRIMMOMATIC=/data/victor/src/Trimmomatic-0.38/trimmomatic-0.38.jar

# New directories
if [ ! -d $TRIMMED ]; then mkdir $TRIMMED; fi
if [ ! -d $TRIMMED/fastqc ]; then mkdir $TRIMMED/fastqc; fi

# Trimming
if [ ! -e $TRIMMED/SRR2029634_reverse_unpaired.fq.gz ]; then
	for i in `seq 27 34`; do
		java -jar $TRIMMOMATIC PE \
			-threads 30 \
			-phred33 \
			$FASTQ/SRR20296$i\_1.fastq \
			$FASTQ/SRR20296$i\_1.fastq \
			$TRIMMED/SRR20296$i\_forward_paired.fq.gz \
			$TRIMMED/SRR20296$i\_forward_unpaired.fq.gz \
			$TRIMMED/SRR20296$i\_reverse_paired.fq.gz \
			$TRIMMED/SRR20296$i\_reverse_unpaired.fq.gz \
			ILLUMINACLIP:$ADAPTERS/TruSeq3-PE-2.fa:2:30:10 \
			ILLUMINACLIP:$ADAPTERS/TruSeq2-PE.fa:2:30:10 \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:40 \
			2>$TRIMMED/SRR20296$i.log 1>$TRIMMED/SRR20296$i.err;
		# Quality control (forward)
		fastqc $TRIMMED/SRR20296$i\_forward_paired.fq.gz \
		       --outdir $TRIMMED/fastqc;
		# Quality control (reverse)
		fastqc $TRIMMED/SRR20296$i\_reverse_paired.fq.gz \
		       --outdir $TRIMMED/fastqc
	done;

fi;

# After trimming, I have checked the fastQC report for both, forward and
# reverse sequences.
# 1 -  It is common to all the reports, to show a bias in the GC
#      composition at the 5' end of the reads. According to fastQC
#      documentation, it looks to be comon to some libraries, but it will not
#      affect the downstream analysis.
# 2 -  Now, a warning appears in the 'Sequence Length Distribution' panel. It
#      was expected because of the trimming. Now, the sequence length variation
#      exists.
# 3 -  The 'Sequence Duplication Levels' warnins persist. It is not worrying
#      because of the low frequence of duplication levels. Less than the 10 %
#      of the reads appears more than once.
# 4 -  There are warnings in all the files in the 'Kmer content' module. It can
#      be dued to the repeated unidentifyied sequences detected.
#
#
# As the quality of the fastq files is acceptable, I will continue with the
# mapping of the files against the Culex quinquefasciatus genome.
