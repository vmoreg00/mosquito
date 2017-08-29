#!/bin/bash
#
#				2017-08-29
#				----------
#
# On 2016-11-22, ipyrad assembled 261 loci with data from at least 2 molestus
# and 2 pipiens individuals. On 2017-05-30, I downloaded reads from other
# mosquito samples, made available by Asgharian et al. (2015), and mapped
# them to the consensus sequences of those loci. However, because most such
# reads came from somewhere else in the genome, the prior probability of any
# read mapping to the loci assembled by ipyrad was very low. As a consequence,
# a large fraction or (partial) mappings were spurious. Then, on 2017-08-04,
# I used the Culex quinquefasciatus reference genome to filter out the reads
# from Asgharian et al. that mapped well anywhere else than in the 261 loci.
# I reduced considerably the number of reads that may map to those loci, and
# that is what I will do now, following the approach in 2017-05-30.

CONSENSUSDIR=../2017-05-30
FASTQDIR=../2017-08-04/fastq
BCFTOOLS=/usr/local/bin/bcftools

# The following are the accession numbers of the raw read from the study by Asgharian et
# al. 2015, Evolutionary genomics of Culex pipiens: global and local adaptations associated
# with climate, life-history traits and anthropogenic factors; Proc. R. Soc. B. 282(1810).
# Here there is some information about the reads, and the newnames that I will use
# to identify each sample:
#
# -----------------------------------------------------------------------------
# Run           Organism        Sampling_site           MegaBases      New name
# -----------------------------------------------------------------------------
# SRR2029627	Culex pipiens	Aleksin Urban A1	2516           PipA1
# SRR2029628	Culex pipiens	Aleksin Suburban A4	4901           PipA4
# SRR2029629	Culex pipiens	Moscow Urban M1 	6618           PipM1
# SRR2029630	Cx. torrentium	Moscow Suburban M2	6535           TorM2
# SRR2029631	Culex pipiens	Sacramento Urban S1	3370           PipS1
# SRR2029632	Culex pipiens	Sacramento Suburban S2	6037           PipS2
# SRR2029633	Culex pipiens	Sacramento Suburban S3	5520           PipS3
# SRR2029634	Cx. torrentium	Moscow Suburban M4	3306           TorM4
# -----------------------------------------------------------------------------
#

SAMPLE=( PipA1 PipA4 PipM1 TorM2 PipS1 PipS2 PipS3 TorM4 )

if [ ! -d bam ]; then mkdir bam; fi

# I realized that I should map the reads in local mode, because many reads
# are expected to overhang the loci selected as references. Using the default
# end-to-end alignment in bowtie2 produces many spurious insertions in the
# edges of the alignments.

for i in ${SAMPLE[@]}; do
   if [ ! -e bam/$i'_sorted.bam' ]; then
      if [ ! -e bam/$i.bam ]; then
         bowtie2 --fast-local \
                 --no-unal \
                 --rg-id $i \
                 --rg SM:$i \
                 --threads 6 \
                 -x $CONSENSUSDIR/consensus \
                 -U $FASTQDIR/$i.fastq 2> bam/$i.log |
         samtools view -Sb - > bam/$i.bam
      fi
      samtools sort bam/$i.bam > bam/$i'_sorted.bam'
      rm bam/$i.bam
      if [ ! -e bam/$i'_sorted.bam.bai' ]; then
         samtools index bam/$i'_sorted.bam'
      fi
   fi
done

if [ ! -e variants.vcf ] && [ ! -e variants.vcf.gz ]; then
   $BCFTOOLS mpileup -f $CONSENSUSDIR/consensus.fasta --read-groups RG.txt --annotate FORMAT/AD -Ou bam/*.bam | \
   $BCFTOOLS call --format-fields GQ,GP --variants-only --multiallelic-caller --prior 0.005 --output variants.vcf
fi

