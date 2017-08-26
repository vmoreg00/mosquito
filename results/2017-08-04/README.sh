#!/bin/bash
#
#				2017-08-04
#				----------
#
# The problem of mapping genomic reads, with multiple origins, on a limited
# subset of loci is the low prior probability of a read mapping to the selected
# loci. This turns into a high false positive rate. In order to increase the
# probability of reads mapping to the selected loci, I should discard the reads
# mapping well elsewhere. I can blast the selected loci (2017-05-30/consensus.fasta)
# to the reference genome of C. quinquefasciatus in order to identify the positions
# of those loci. Then, I will map the reads, with global settings, to that reference
# genome, and select the ones that map around the loci identified, together with the
# ones that did not map anywhere, just in case.

LOCI=../2017-05-30/consensus.fasta
READSDIR=../2017-05-30/trimmed

if [ ! -e Culex.fa.nhr ]; then
   if [ ! -e Culex.fa ]; then
      if [ ! -e Culex_quinquefasciatus.CpipJ2.dna_sm.toplevel.fa.gz ]; then
         wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-36/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.dna_sm.toplevel.fa.gz
      fi
      gunzip Culex_quinquefasciatus.CpipJ2.dna_sm.toplevel.fa.gz
      mv Culex_quinquefasciatus.CpipJ2.dna_sm.toplevel.fa Culex.fa
   fi
   makeblastdb -in Culex.fa -dbtype nucl -logfile makeblastdb.log
   # rm Culex.fa
fi

if [ ! -e loci.blast ]; then
   # I tried with megablast, dc-megablast and blastn. The numbers of queries with
   # at least one hit were 144, 161, and 260, respectively, out of 261 queries.
   blastn -query $LOCI \
          -db Culex.fa \
          -task blastn \
          -out loci.blast \
          -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send sstrand length pident score bitscore evalue'
fi

# No more than 115 out of 260 loci with at least one hit align through more than 90%
# of their length somewhere in the genome. Loci that do not find their place in the
# genome may still be valid, though. The loci with more than one hit are not, and they
# must be excluded.
#
# The goal is to identify the regions where valid loci exist. Only reads mapping to those
# loci, or not mapping anywhere else in the reference genome will be used to genotype
# the loci.
#
# Valid loci have unique, almost full length alignments to non-mitochondrial contigs. The
# uniqueness can be assessed by the difference (or ratio) between the e-values (or scores)
# of the first and the second hits. I consider unique a best hit with a score at least 1.5
# times larger than the second best.

if [ ! -e valid_regions.bed ]; then
   # I need to sort the blast hits first by qseqid, and then by bitscore in reverse order.
   LC_ALL=C sort -k 1.7n,1 -k 13gr,13 -k 3n,3 -k 5.13n,5 -k 7n,7 loci.blast | \
   gawk '{
      HIT[$1]++
      BITSCORE[$1, HIT[$1]] = $13
      QALIGNED[$1, HIT[$1]] = $4 - $3 + 1
      QLEN[$1] = $2
      SSEQID[$1, HIT[$1]] = $5
      if ($9 == "plus") {
         SSTART[$1, HIT[$1]] = $7
         SEND[$1, HIT[$1]] = $8
      } else {
         SSTART[$1, HIT[$1]] = $8
         SEND[$1, HIT[$1]] = $7
      }
   }END{
      for (QSEQID in HIT) {
         if (HIT[QSEQID] > 1) {
            if ((BITSCORE[QSEQID, 1] / BITSCORE[QSEQID, 2] > 1.5) && (QALIGNED[QSEQID, 1] > 0.8 * QLEN[QSEQID]) && (SSEQID[QSEQID, 1] !~ "Mt")) {
               print SSEQID[QSEQID, 1] "\t" SSTART[QSEQID, 1] - 1 "\t" SEND[QSEQID, 1] "\t" QSEQID
            }
         } else {
            if ((QALIGNED[QSEQID, 1] > 0.8 * QLEN[QSEQID]) && (SSEQID[QSEQID, 1] !~ "Mt")) {
               print SSEQID[QSEQID, 1] "\t" SSTART[QSEQID, 1] - 1 "\t" SEND[QSEQID, 1] "\t" QSEQID
            }
         }
      }
   }' | LC_ALL=C sort -k 1.12g,1 -k 2n,2 -k 3n,3 > valid_regions.bed
fi

# I notice some overlaps of valid regions that should correspond to different loci.
# locus_555106 locus_367966 map to opposite strands of the same region. When I align
# the consensus sequences from the two loci with exonerate, I obtain this alignment:
#
# C4 Alignment:
# ------------
#          Query: locus_367966 [revcomp]
#         Target: locus_555106
#          Model: ungapped:dna2dna
#      Raw score: 561
#    Query range: 247 -> 3
#   Target range: 134 -> 378
#
#  247 : TGATCTTTCTTTGGCCCAACTCCCCTTCCCTCTCTAATTTCAATTGTCATTGCTCGACTTACTA : 184
#        ||||| ||||  || |||   || |  ||:|| |  ||:|| || ||||| |  ||    || |
#  135 : TGATCCTTCTCAGGKCCACAACCGCGACCYTCGCGGATYTCGATRGTCATCGAACGGGACACSA : 198
#
#  183 : CATCACGAGAAGCCAAATCTTTTGCCTTTGGCGCATAACGTTCCATAAACTTCTCACCCTGAGA : 120
#        | ||||| || ||:||||| || ||:   || ||||| || ||||| ||    |||||||||||
#  199 : CGTCACGGGAGGCYAAATCCTTCGCYACCGGTGCATATCGCTCCATGaaacgttcaccctgaga : 262
#
#  119 : ATTAACGAGGTACCCACCTTCACCGCGGCATCCTTCCGTCATCAAGCACCCTGAGCCATATATT :  56
#         || |  ||||| |||||||| || |||||||| || || ||||  || || |  || || |||
#  263 : gttgatcaggtagccaccttctcctcggcatccctcggtgatcagacaaccggcaccgtagatt : 326
#
#   55 : CCTGTTGGATGAAATTGCACAAAYTCCATATCTTCAAGTGGCAACCCAGCTC :   4
#        || || || || || ||||||||:|||| ||| ||    ||||  ||||| |
#  327 : cccgtcgggtggaactgcacaaactccagatcctccgagggcagtccagcac : 378
#
# vulgar: locus_367966 247 3 - locus_555106 134 378 + 561 M 244 244
#
# Other apparent overlaps are just due to the two loci being adjacent to each other.
# In principle, in this step I only need to determine to what regions of a reference
# genome reads should map in order to believe them as evidence of the genotypes at
# certain loci. I do not need to exclude loci now. Reads that seem to be generated
# by a valid region may map to both loci, and get discarded then.
#
# Now, let's map the available reads to the reference genome and extract the reads
# that map either to the valid regions or nowhere in the genome. A strict mapping
# setting is necessary to exclude from further analysis only the reads that we can
# confidently attribute to loci that are not of interest. End-to-end alignment is
# in order.

if [ ! -e Culex.1.bt2 ]; then
   bowtie2-build Culex.fa Culex
fi

if [ ! -d bam ]; then mkdir bam; fi

# Recall the origin of the reads, published by Asgharian et al. (2015):
#
# ---------------------------------------------------------------
# Run           Organism        Sampling_site           MegaBases
# ---------------------------------------------------------------
# SRR2029627    Culex pipiens   Aleksin Urban A1        2516
# SRR2029628    Culex pipiens   Aleksin Suburban A4     4901
# SRR2029629    Culex pipiens   Moscow Urban M1         6618
# SRR2029630    Cx. torrentium  Moscow Suburban M2      6535
# SRR2029631    Culex pipiens   Sacramento Urban S1     3370
# SRR2029632    Culex pipiens   Sacramento Suburban S2  6037
# SRR2029633    Culex pipiens   Sacramento Suburban S3  5520
# SRR2029634    Cx. torrentium  Moscow Suburban M4      3306
# ---------------------------------------------------------------
#
# Following the analysis in 2017-07-19, I will change the original group names:

ORIGINAL=( SRR2029627 SRR2029628 SRR2029629 SRR2029630
           SRR2029631 SRR2029632 SRR2029633 SRR2029634 )

NEWNAME=( PipA1 PipA4 PipM1 TorM2 PipS1  PipS2  PipS3 TorM4 )


for i in 0 1 2 3 4 5 6 7; do
   if [ ! -e bam/${NEWNAME[$i]}'_sorted.bam' ]; then
      if [ ! -e bam/${NEWNAME[$i]}.bam ]; then
         bowtie2 --sensitive \
                 --rg-id ${NEWNAME[$i]} \
                 --rg SM:${NEWNAME[$i]} \
                 --threads 6 \
                 -x Culex \
                 -U $READSDIR/${ORIGINAL[$i]}.forward.fastq,$READSDIR/${ORIGINAL[$i]}.reverse.fastq,$READSDIR/${ORIGINAL[$i]}.merged.fastq 2> bam/${NEWNAME[$i]}.log |
         samtools view -Sb - > bam/${NEWNAME[$i]}.bam
      fi
      samtools sort bam/${NEWNAME[$i]}.bam > bam/${NEWNAME[$i]}'_sorted.bam'
      rm bam/${NEWNAME[$i]}.bam
      if [ ! -e bam/${NEWNAME[$i]}'_sorted.bam.bai' ]; then
         samtools index bam/${NEWNAME[$i]}'_sorted.bam'
      fi
   fi
done

if [ ! -e bam/summary.txt ]; then
   echo -e "Sample\t        Unaligned\t     Aligned_once\tAligned_ambiguous\tAlignment_rate" > bam/summary.txt
   for i in 4 5 6 0 1 2 3 7; do
      gawk -v SAMPLE=${NEWNAME[$i]} '(/aligned 0 times/){
         UNALIGNED = $1 " " $2
      }(/aligned exactly 1 time/){
         ONCE = $1 " " $2
      }(/aligned >1 times/){
         MULTIPLE = $1 " " $2
      }(/overall alignment rate/){
         RATE = $1
      }END{
         printf("%s\t% 17s\t% 17s\t% 17s\t%s\n", SAMPLE, UNALIGNED, ONCE, MULTIPLE, RATE)
      }' bam/${NEWNAME[$i]}.log >> bam/summary.txt
   done
fi

# ----------------------------------------------------------------------------------------------------
# Sample	        Unaligned	     Aligned_once	Aligned_ambiguous	Alignment_rate
# ----------------------------------------------------------------------------------------------------
# PipS1 	 8981600 (30.95%)	 8338628 (28.74%)	11696049 (40.31%)	        69.05%
# PipS2 	11377823 (25.13%)	16065682 (35.48%)	17840750 (39.40%)	        74.87%
# PipS3 	13288590 (32.36%)	11939942 (29.07%)	15839026 (38.57%)	        67.64%
# PipA1 	 6552769 (34.45%)	 6503921 (34.20%)	 5962175 (31.35%)	        65.55%
# PipA4 	13982465 (37.30%)	 7879847 (21.02%)	15619511 (41.67%)	        62.70%
# PipM1 	18846726 (35.33%)	14952687 (28.03%)	19543652 (36.64%)	        64.67%
# TorM2 	28411160 (54.67%)	 9463321 (18.21%)	14096170 (27.12%)	        45.33%
# TorM4 	16202316 (57.99%)	 3722156 (13.32%)	 8015475 (28.69%)	        42.01%
# ----------------------------------------------------------------------------------------------------
#
# Now, I should select the reads that are either unaligned or aligned confidently
# to the valid loci.

if [ ! -d fastq ]; then mkdir fastq; fi

for i in 0 1 2 3 4 5 6 7; do
   if [ ! -e fastq/${NEWNAME[$i]}.fastq ]; then
      samtools view -f 4 bam/${NEWNAME[$i]}'_sorted.bam' | \
      gawk -v DESCRIPTION="not_placed" -f sam2fastq.awk > fastq/${NEWNAME[$i]}.fastq

      samtools view -q 30 -L valid_regions.bed bam/${NEWNAME[$i]}'_sorted.bam' | \
      gawk -v DESCRIPTION="from_valid" -f sam2fastq.awk >> fastq/${NEWNAME[$i]}.fastq
   fi
done
