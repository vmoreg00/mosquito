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

if [ ! -d mappeo ]; then mkdir mappeo; fi

if [ ! -e culex.1.bt2 ]; then
   bowtie2-build reference.fa culex
fi

for i in `seq 0 16`; do
   if [ ! -e mappeo/${LISTA[$i]}'_r1'.sam ]; then
      bowtie2 --local \
              --very-sensitive \
              -x culex \
              -U $FASTQDIR/${LISTA[$i]}'_R1_trimmed.fastq' \
              -S mappeo/${LISTA[$i]}'_r1'.sam \
              --rg-id ${LISTA[$i]}
   fi
   if [ ! -e mappeo/${LISTA[$i]}'_r2'.sam ]; then
      bowtie2 --local \
              --very-sensitive \
              -x culex \
              -U $FASTQDIR/${LISTA[$i]}'_R2_trimmed.fastq' \
              -S mappeo/${LISTA[$i]}'r2'.sam \
              --rg-id ${LISTA[$i]}
   fi
   if [ ! -e mappeo/${LISTA[$i]}.sam ]; then
      bowtie2 --local \
              --very-sensitive \
              -x culex \
              -U $FASTQDIR/${LISTA[$i]}'_trimmed.fastq' \
              -S mappeo/${LISTA[$i]}.sam \
              --rg-id ${LISTA[$i]}
   fi

done

if [ ! -d BAM ]; then mkdir BAM; fi

for i in `seq 0 16`; do
	if [ ! -e BAM/${LISTA[$i]}'_r1.bam' ]; then
		samtools view -bS mappeo/${LISTA[$i]}'_r1'.sam > BAM/${LISTA[$i]}'_r1.bam'
	fi
	if [ ! -e BAM/${LISTA[$i]}'_r2.bam' ]; then
                samtools view -bS mappeo/${LISTA[$i]}'_r2'.sam > BAM/${LISTA[$i]}'_r2.bam'
	fi
        if [ ! -e BAM/${LISTA[$i]}'_merged.bam' ]; then
                samtools view -bS mappeo/${LISTA[$i]}.sam > BAM/${LISTA[$i]}'_merged.bam'
        fi
done

