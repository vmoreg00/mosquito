#!/bin/bash
#
#				2017-05-30
#				----------
#
# In folder 2016-11-22 ipyrad assembled 261 loci with at least 2 molestus
# and 2 pipiens individuals. Here, I want to build a reference genome with
# the consenus sequences from each locus. Then, I will download available
# reads and map them to this reference to increase the dataset.
#

ALLELES=../2016-11-22/culex/min2_outfiles/min2.alleles.loci
SRATOOL=/home/joiglu/bin/sratoolkit.2.8.2-1-ubuntu64/bin

if [ ! -d loci ]; then mkdir loci; fi
if [ `ls -1 loci | wc -l` -eq 0 ]; then
   gawk '(/^[^\/]/){
      SEQ[$1]  = $2
   }(/^\/\//){
      split($NF,A,/\|/)
      FILE = "loci/locus_" A[2] ".fasta"
      for (name in SEQ) {
         print ">" name  >> FILE
         print SEQ[name] >> FILE
      }
      close(FILE)
      delete(SEQ)
      delete(A)
   }' $ALLELES
fi

if [ ! -e consensus.fasta ]; then
   for locus in `ls -1 loci/`; do
#      LOCUS=`basename ${locus:6} .fasta`
      consambig -sequence loci/$locus -outseq z1.fasta -name `basename $locus .fasta`
      cat z1.fasta >> consensus.fasta
      rm z1.fasta
   done
fi

if [ ! -e consensus.1.bt2 ]; then
   bowtie2-build consensus.fasta consensus
fi

# The following are the accession numbers of the raw read from the study by Asgharian et
# al. 2015, Evolutionary genomics of Culex pipiens: global and local adaptations associated
# with climate, life-history traits and anthropogenic factors; Proc. R. Soc. B. 282(1810).
# Here there is some information about the reads:
#
# ---------------------------------------------------------------
# Run           Organism        Sampling_site           MegaBases
# ---------------------------------------------------------------
# SRR2029627	Culex pipiens	Aleksin Urban A1	2516
# SRR2029628	Culex pipiens	Aleksin Suburban A4	4901
# SRR2029629	Culex pipiens	Moscow Urban M1 	6618
# SRR2029630	Cx. torrentium	Moscow Suburban M2	6535
# SRR2029631	Culex pipiens	Sacramento Urban S1	3370
# SRR2029632	Culex pipiens	Sacramento Suburban S2	6037
# SRR2029633	Culex pipiens	Sacramento Suburban S3	5520
# SRR2029634	Cx. torrentium	Moscow Suburban M4	3306
# ---------------------------------------------------------------
#
# There is hardly any information on how the authors processed the reads. It may be necessary
# to trim and filter the reads before mapping them.

if [ ! -d fastq ]; then mkdir fastq; fi
if [ ! -d fastqc ]; then mkdir fastqc; fi

for i in SRR2029627 SRR2029628 SRR2029629 SRR2029630 \
         SRR2029631 SRR2029632 SRR2029633 SRR2029634; do
   if [ ! -e fastq/$i'_1.fastq' ]; then
      $SRATOOL/fastq-dump -F --split-files --outdir ./fastq $i
   fi
   if [ ! -e fastqc/$i'_1_fastqc.html' ]; then
      fastqc fastq/$i'_1.fastq' --outdir ./fastqc &
   fi
   if [ ! -e fastqc/$i'_2_fastqc.html' ]; then
      fastqc fastq/$i'_2.fastq' --outdir ./fastqc &
   fi
done

wait

# Some fastq files have overrepresented sequences, which always match some adapter or
# contaminant sequence. It is necessary to inspect manually what adapters those sequences
# belong to. For example, the hits to TruSeq Adapter with Index 7 and those to Illumina
# Multiplexing PCR Primer 2.01 in forward reads actually correspond to the same sequence,
# which seems to be the TruSeq Adapter with index 7. This adapter appears in forward reads
# in the  same sense in which it is defined in the Illumina documentation. Note that the
# overrepresented sequences identified as possibly coming from Illumina's Multiplexing PCR
# Primer 2.01 would match the reverse complementary of that primer, but aligns directly
# to the TruSeq Adapter.
#
# On reverse reads, the overrepresented sequences are identified as possibly coming from
# Illumina Single End PCR Primer 1, although they would match it in the opposite strand.
# In fact, the overrepresented string can be extended to the full length of the reverse
# complement of the Universal TruSeq Adapter.
#
# Thus, the library must have been something like this:
#
#          Univ.TruSeq              Genomic            TruSeq In. X
#    5'|-------------------|---------------------|---------------------|3`
#
# The index corresponding to each sample is the following:
#
# -------------------------------------------------------------------------------------
# Sample    	 Index	Sequence
# -------------------------------------------------------------------------------------
# SRR2029627	     7	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
# SRR2029628	    10  GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
# SRR2029629	     1	GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
# SRR2029630	     6	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
# SRR2029631	     3	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
# SRR2029632	     4	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
# SRR2029633	     5	GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
# SRR2029634	     2	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
# -------------------------------------------------------------------------------------
#
# The sequence of the common adapter that appears in reverse reads is this:
#
#    AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
#
# But what gets sequenced is its reverse complementary, since reverse reads are reading
# the bottom strand. So it should be searched as this:
#
#    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#
# However, I will search that adapter on the reverse reads that did not get assembled
# with PEAR (I don't expect any in the assembled ones). And recall that PEAR reverse-
# complements the reverse reads. Therefore, I will search the top strand of the adapter
# and in 5-prime, instead of 3-prime.

if [ ! -d merged ]; then mkdir merged; fi
for i in `seq 27 34`; do
   if [ ! -e merged/SRR20296$i.log ]; then
      pear -f fastq/SRR20296$i'_1.fastq' \
           -r fastq/SRR20296$i'_2.fastq' \
           -o merged/SRR20296$i > merged/SRR20296$i.log
   fi
done

if [ ! -d trimmed ]; then mkdir trimmed; fi
for i in `seq 27 34`; do
   if [ ! -e trimmed/SRR20296$i.forward.log ]; then
      cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               --max-n=5 \
               --minimum-length=50 \
               -o trimmed/SRR20296$i.forward.fastq \
               merged/SRR20296$i.unassembled.forward.fastq > trimmed/SRR20296$i.forward.log &
   fi
   if [ ! -e trimmed/SRR20296$i.reverse.log ]; then
      cutadapt -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
               --max-n=5 \
               --minimum-length=50 \
               --error-rate=0.2 \
               -o trimmed/SRR20296$i.reverse.fastq \
               merged/SRR20296$i.unassembled.reverse.fastq > trimmed/SRR20296$i.reverse.log &
   fi
   if [ ! -e trimmed/SRR20296$i.merged.log ]; then
      cutadapt -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
               -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               --max-n=5 \
               --minimum-length=50 \
               -o trimmed/SRR20296$i.merged.fastq \
               merged/SRR20296$i.assembled.fastq > trimmed/SRR20296$i.merged.log &
   fi
   wait
done

if [ ! -d bam ]; then mkdir bam; fi

# I realized that I should map the reads in local mode, because many reads
# are expected to overhang the loci selected as references. Using the default
# end-to-end alignment in bowtie2 produces many spurious insertions in the
# edges of the alignments.

for i in `seq 27 34`; do
   if [ ! -e bam/SRR20296$i'_sorted.bam' ]; then
      if [ ! -e bam/SRR20296$i.bam ]; then
         bowtie2 --fast-local \
                 --no-unal \
                 --rg-id SRR20296$i \
                 --rg SM:SRR20296$i \
                 --threads 6 \
                 -x consensus \
                 -U trimmed/SRR20296$i.forward.fastq,trimmed/SRR20296$i.reverse.fastq,trimmed/SRR20296$i.merged.fastq 2> bam/SRR20296$i.log |
         samtools view -Sb - > bam/SRR20296$i.bam
      fi
      samtools sort bam/SRR20296$i.bam > bam/SRR20296$i'_sorted.bam'
      rm bam/SRR20296$i.bam
      if [ ! -e bam/SRR20296$i'_sorted.bam.bai' ]; then
         samtools index bam/SRR20296$i'_sorted.bam'
      fi
   fi
done
