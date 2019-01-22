# Test of correlation between the qualities of freebayes and GATK
# date: 2018-11-29
# author: Victor Moreno Gonzalez

# input files
fb=/data/victor/mosquito/results/2018-07-25/freebayes_test/vcf/culex_merge.vcf
gatk=/data/victor/mosquito/results/2018-07-25/gatk_test/vcf/culex_cohort.vcf

# extract the positions and qualities of each vcf
grep -v "^#" $gatk | cut -f 1,2,6,5 | awk '{print "G", $0}' > gatk_qual.csv
grep -v "^#" $fb | cut -f 1,2,6,5 | awk '{print "F", $0}' > fb_qual.csv
cat gatk_qual.csv fb_qual.csv | sort -k2r,2 -k3 > both_qual.sorted.csv
