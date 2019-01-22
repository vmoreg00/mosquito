# script to launch crisp while freebayes finishes. This command lines
# should go into the README file of the parent directory.
REFGENOME=/data/victor/mosquito/data/refgenome
BAM=/data/victor/mosquito/data/asgharian/bam/bwa
VCF=/data/victor/mosquito/data/asgharian/vcf

sample=(1.molA1 2.pipA4 3.molM1 4.torM2 5.molS1 6.mixS2 7.mixS3 8.torM4) # Names
ploidy=(41 41 52 56 30 26 41 41); # Genome copies in each sample

if [ ! -d $VCF/crisp ]; then
        mkdir $VCF/crisp;
fi;

if [ ! -e bam.list ]; then
	for i in `seq 0 7`; do
        	echo "$BAM/${sample[$i]}.sort.RG.markdup.bam PS=${ploidy[$i]}" \
		     >> bam.list
	done;
fi;

# variant calling is performed over the merged bam made during the freebayes
# variant calling.
# Constraints (minc, mbq, mmq, filterreads) have been modifyied to detect
# some more variants. With the default parameters, the number of SNVs were too
# low.
if [ ! -e $VCF/crisp/culex_cohort.vcf ]; then
	crisp --bams bam.list \
	      --ref $REFGENOME/CulQui.fna \
	      --VCF $VCF/crisp/culex_cohort.vcf \
	      --minc 1 \
	      --mbq 5 \
	      --mmq 10 \
	      --perms 30000 \
	      --filterreads 0 \
	      --qvoffset 33 \
	      --EM 1 \
	      --verbose 1 > $VCF/crisp/culex_cohort.log
fi;
