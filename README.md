# Mosquito

This is a re-analysis of some preliminary sequence data from Culex pipiens
form pipiens and form molestus. Individual samples were sequenced at a very
shallow depth. Previous analyses showed a quite low level of mapping success
to the reference genome of Culex quinquefasciatus. Chimeric reads are suspected
to be frequent.

The goal is to compare the two ecotypes.

#2016-07-26
We subsampled molestus reads and used pyrad to cluster the reads under different
similarity thresholds. The results suggest that a similarity threshold of 95%
provides the highest coverage per cluster, before dropping abruptly. However, the
fact that the number of clusters increase steadily suggests extensive paralogy
problems.

#2016-07-25
Here we generate an alternative vcf file with Freebayes, instead of samtools'
mpileup (2016-07-12). Freebayes is much more flexible, newer, and probably better.

#2016-07-21
First attempt to run pyrad. Very few reads, available, though. In any case, this is
deprecated, since we should be using ipyrad.

#2016-07-12
We used samtools mpileup to generate the vcf. A first look at the distribution of
configurations suggests that the low coverage among pipiens samples precludes any
attempt at guessing genotypes, and it seems that in most sites only one or very few
individuals have any coverage. The situation is better for the molestus samples.
I need to revise the version of samtools used, which may affect the vcf.

#2016-07-08
We used bowtie2 to map the merged reads to the reference genome of Culex quinquefasciatus.
They map better under local settings. Around 94% of molestus reads are mapped (only
25% unambiguously, though), while only 64% of pipiens reads map at all (only 16%
unambiguously). A possibility is to map the original paired ends as independent
reads, to allow for chimeric reads to map better. This would involve changing the
names of the reads, and it is not done, yet.

#2016-07-06b
Now we check for the presence of adapters and trim them with Cutadapt. Read pairs
that have been successfully merged by PEAR are not expected to contain adapters
in 3' ends. Indeed, no more than 0.1% of reads seem to have (spurious?) traces of
adapters.

#2016-07-06
We run SABRE to demultiplex the original fastq files, and PEAR, to merge the
paired reads that overlap. Both molestus and pipiens samples are processed. Less
than 10% of reads fail to be assinged to any sample. More than 95% of identified
pairs of reads merge. It is likely that many or most of those pairs that do not
merge fail to merge because of lower quality. It is safe to disregard non-merged
reads and work only with the assembled subset from now on.

