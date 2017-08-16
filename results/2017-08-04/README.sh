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

if [ ! -e valid_regions.txt ]; then
   # I need to sort the blast hits first by qseqid, and then by bitscore in reverse order.
   sort -k 1.7n,1 -k 13gr,13 -k 3n,3 -k 5.13n,5 -k 7n,7 loci.blast | \
   gawk '{
      HIT[$1]++
      BITSCORE[$1, HIT[$1]] = $13
      QALIGNED[$1, HIT[$1]] = $4 - $3 + 1
      QLEN[$1] = $2
      SSEQID[$1, HIT[$1]] = $5
      if ($9 ~ plus) {
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
               print QSEQID "\t" SSEQID[QSEQID, 1] "\t" SSTART[QSEQID, 1] - 200 "\t" SEND[QSEQID, 1] + 200
            }
         } else {
            if ((QALIGNED[QSEQID, 1] > 0.8 * QLEN[QSEQID]) && (SSEQID[QSEQID, 1] !~ "Mt")) {
               print QSEQID "\t" SSEQID[QSEQID, 1] "\t" SSTART[QSEQID, 1] - 200 "\t" SEND[QSEQID, 1] + 200
            }
         }
      }
   }' > valid_regions.txt
fi
