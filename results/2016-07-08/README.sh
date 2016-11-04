#!/bin/bash
#
#				2016-07-08
#				----------
#
# Mapping reads with bowtie2 to Culex quinquefasciatus genome. The reads to map are
# the fastq files that had been merged and cleaned with cutadapt. They should be in
# ../2016-07-06b.

FASTQDIR='../2016-07-06b'

if [ ! -e reference.fa ]; then
	wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz
	mv Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz reference.fa.gz
	gunzip reference.fa.gz
fi

LISTA=(PipFe1 PipFe2 PipFe3 PipFe4 PipFe5 PipFe6 PipMa1 PipMa2 \
       PipMa3 PipMa4 PipMa5 PipMa6 Mol01 Mol02 Mol03 Mol04 Mol05)

if [ ! -e culex.1.bt2 ]; then
   bowtie2-build reference.fa culex
fi

# The first time we made the mapping wasn't very successful so in this time we decided to make a local alignment.
# Also we  processed the unassembled reads like if they were pair-ends but the mapping was not good again, so
# finally we decided to process the reads, assembled and unassembled, as independent ones. The rationale
# is that if reads map better with local alignment than with a global one (presumably due to
# many reads being chimeric), they will also map better as single ends than as paired ends, even
# if they are paired in nature.
#
# This introduces the problem of the names. Paired reads have the same name in the fastq files, only
# distinguished by the second word of the name line in the fastq format. Only the first word is retained
# in the sam output. Thus, two paired reads would appear in sam as having two primary alignments, if
# bowtie2 allows for two reads with the same name at all.
#
# The natural thing to do is to change the names of the second reads. One option is to use bbreformat.sh,
# from the bbmap package. This program allows to add ' /1' and ' /2' to the names of the reads. But that
# does not help, because of the white space. Another option is to substitute spaces by underscores in the
# names.
#
# It is desirable to run in parallele not only the mapping but also the re-writing of the fastq files
# and the post-mapping creation of the bam file. Thus, it is best to outsource the processing of each
# sample to a different script (it could have been a function), that we can call from here.
#
# Note that the execution of the auxiliary script is made conditional on the non-existence of the final
# bam output, so that intermediate files, such as the created fastq files and the sam files can be deleted.

PROC=`grep -P '^processor' /proc/cpuinfo | wc -l`
if [ $PROC -gt 17 ]; then
   for i in `seq 0 16`; do
      if [ ! -e ${LISTA[$i]}'.bam' ]; then
#        ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} &
         ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} CleanUp &

      fi
   done
   wait
else
   for i in 0 1 2 3 4 5; do
      if [ ! -e ${LISTA[$i]}'.bam' ]; then
#        ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} &
         ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} CleanUp &
      fi
   done
   wait
   for i in 6 7 8 9 10 11; do
      if [ ! -e ${LISTA[$i]}'.bam' ]; then
#        ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} &
         ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} CleanUp &

      fi
   done
   wait
   for i in 12 13 14 15 16; do
      if [ ! -e ${LISTA[$i]}'.bam' ]; then
#        ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} &
         ./map_and_make_bam.sh $FASTQDIR ${LISTA[$i]} CleanUp &

      fi
   done
   wait
fi
# After everything is done, we run the summary:
if [ ! -e summary_map.txt ]; then
   if [ ! -e archivos.txt ]; then
      ls -1 *.maplog > archivos.txt
   fi
   ./summary.py archivos.txt
    rm archivos.txt
fi

# Actually, the merged reads are the ones that got mapped. We did not change the names of
# individual ends, nor map them independently. Overall, the alignment rate seems good, and
# specially so for the molestus samples, presumably less affected by chimeras. However, the
# portion of reads mapped unambiguously is still very low.
