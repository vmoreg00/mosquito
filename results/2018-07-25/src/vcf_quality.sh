##############################################################################
#                                                                            #
# vcf_quality.sh                                                             #
#                                                                            #
# Script to check the quality of the variant calling performed by both,      #
# GATK and freebayes                                                         #
#                                                                            #
# Author: Victor Moreno-Gonzalez                                             #
# Date: 09-10-2018                                                           #
#                                                                            #
##############################################################################

VCF=/data/victor/mosquito/data/asgharian/vcf
REFGENOME=/data/victor/mosquito/data/refgenome
OUTDIR=/data/victor/mosquito/results/2018-07-25/stats

SNP=/data/victor/mosquito/results/2018-07-25/src/snp_qual_dist.py

if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR;
fi;

# GATK variant calling quality
if [ ! -d $OUTDIR/GATK ]; then
	mkdir $OUTDIR/GATK;
fi;

for i in `seq 27 34`; do
	if [ ! -e $OUTDIR/GATK/SRR20296$i.stats ]; then
		bcftools stats -F $REFGENOME/CulQui.fna \
			$VCF/GATK/SRR20296$i.vcf > $OUTDIR/GATK/SRR20296$i.stats &
	fi;
	if [ ! -e $OUTDIR/GATK/SRR20296$i\_SNPs.png ]; then
		$SNP $VCF/GATK/SRR20296$i.vcf $OUTDIR/GATK/SRR20296$i\_SNPs &
	fi;
done;
wait

for i in `seq 27 34`; do
	if [ ! -e $OUTDIR/GATK/SRR20296$i.pdf ]; then
		plot-vcfstats --prefix $OUTDIR/GATK/tmp/ \
			--title SRR20296$i --main-title 20296$i \
			$OUTDIR/GATK/SRR20296$i.stats;
		mv $OUTDIR/GATK/tmp/summary.pdf $OUTDIR/GATK/SRR20296$i.pdf;
	fi;
done;

if [ -d $OUTDIR/GATK/tmp ]; then
	rm -r $OUTDIR/GATK/tmp;
fi;


# freebayes variant calling quality
if [ ! -d $OUTDIR/FB ]; then
	mkdir $OUTDIR/FB;
fi;

for i in `seq 27 34`; do
        if [ ! -e $OUTDIR/FB/SRR20296$i.stats ]; then
                bcftools stats -F $REFGENOME/CulQui.fna \
                        $VCF/FB/SRR20296$i.vcf > $OUTDIR/FB/SRR20296$i.stats &
        fi;
        if [ ! -e $OUTDIR/FB/SRR20296$i\_SNPs.png ]; then
                $SNP $VCF/FB/SRR20296$i.vcf $OUTDIR/FB/SRR20296$i\_SNPs &
        fi;
done;
wait

for i in `seq 27 34`; do
        if [ ! -e $OUTDIR/FB/SRR20296$i.pdf ]; then
                plot-vcfstats --prefix $OUTDIR/FB/tmp/ \
                        --title SRR20296$i --main-title 20296$i \
                        $OUTDIR/FB/SRR20296$i.stats;
                mv $OUTDIR/FB/tmp/summary.pdf $OUTDIR/FB/SRR20296$i.pdf;
        fi;
done;

if [ -d $OUTDIR/FB/tmp ]; then
        rm -r $OUTDIR/FB/tmp;
fi;

# filtered variant calling quality
if [ ! -d $OUTDIR/filt ]; then
        mkdir $OUTDIR/filt;
fi;

for i in `seq 28 34`; do
        if [ ! -e $OUTDIR/filt/SRR20296$i.stats ]; then
                bcftools stats -F $REFGENOME/CulQui.fna \
                        $VCF/filt/SRR20296$i.vcf > $OUTDIR/filt/SRR20296$i.stats &
        fi;
        if [ ! -e $OUTDIR/filt/SRR20296$i\_SNPs.png ]; then
                $SNP $VCF/filt/SRR20296$i.vcf $OUTDIR/filt/SRR20296$i\_SNPs &
        fi;
done;
wait

for i in `seq 28 34`; do
        if [ ! -e $OUTDIR/filt/SRR20296$i.pdf ]; then
                plot-vcfstats --prefix $OUTDIR/filt/tmp/ \
                        --title SRR20296$i --main-title 20296$i \
                        $OUTDIR/filt/SRR20296$i.stats;
                mv $OUTDIR/filt/tmp/summary.pdf $OUTDIR/filt/SRR20296$i.pdf;
        fi;
done;

if [ -d $OUTDIR/filt/tmp ]; then
        rm -r $OUTDIR/filt/tmp;
fi;

