# Mosquito

This is a re-analysis of some preliminary sequence data from Culex pipiens
form pipiens and form molestus. Individual samples were sequenced at a very
shallow depth. Previous analyses showed a quite low level of mapping success
to the reference genome of Culex quinquefasciatus. Chimeric reads are suspected
to be frequent.

The goal is to compare the two ecotypes.

#2016-07-26
We decided to run pyrad with different values of Wclust: 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99. Also, we made a table to be able to visualize the results.

#2016-07-25
Another way to evaluate mapping is  using a programme called Freebayes. It is a a bayesian genetic variant detector to find SNPs, indels, MNPs and complex events. It has a algorithm better than, for example samtools, because it calls variants based on the literal sequences of reads aligned to a particular target, not their precise alignment.

#2016-07-21
Run pyrad  with only the assembled reads.

#2016-07-12
To be able to visualize SNPs and indels we converted BAM files to VCF. First of all, we downloaded and indexed the reference. After that, the samples were sorted and indexed. Once you have it, you will be able to convert to VCF and to analyze any genotypic variation.

#2016-07-08
We used bowtie2 to map the reads to the reference genome of Culex quinquefasciatus. By mistake, we did not use the trimmed reads, but

We wanted to check why our mapping wasn't good enough. We used Bowtie2 to map our 'independent'-reads against the reference. Furthermore, we converted sam to bam.

#2016-07-06b
Now we check for the presence of adapters and trim them with Cutadapt. Read pairs
that have been successfully merged by PEAR are not expected to contain adapters
in 3' ends.

#2016-07-06
We run SABRE to demultiplex the original fastq files, and PEAR, to merge the
paired reads that overlap. Both molestus and pipiens samples are processed.




