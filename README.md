# Mosquito

This is a re-analysis of some preliminary sequence data from Culex pipiens
form pipiens and form molestus. Individual samples were sequenced at a very
shallow depth. Previous analyses showed a quite low level of mapping success
to the reference genome of Culex quinquefasciatus. Chimeric reads are suspected
to be frequent.

The goal is to compare the two ecotypes.

#2016-07-06
We run SABRE to demultiplex the original fastq files, and PEAR, to merge the
paired reads that overlap. Both molestus and pipiens samples are processed.

#2016-07-06b
Now we check for the presence of adapters and trim them with Cutadapt. Read pairs
that have been successfully merged by PEAR are not expected to contain adapters
in 3' ends.

#2016-07-08
We used bowtie2 to map the reads to the reference genome of Culex quinquefasciatus.
By mistake, we did not use the trimmed reads, but

We wanted to check why our mapping wasn't good enough. We used Bowtie2 to map our
'independent'-reads against the reference. Furthermore, we converted sam to bam.

#2016-07-12
To be able to visualize SNPs and indels we converted BAM files to VCF. First of all, we downloaded and indexed the reference. After that, the samples were sorted and indexed. Once you have it, wou will be able to convert to VCF and to analize any genotypic variation
