#!/bin/bash
#$ -N ipyrad_culex
#$ -m ea
#$ -M joiglu@uv.es
#$ -e readme.err
#$ -o readme.log
#$ -wd /scratch/jo/joiglu/2016-11-22
#$ -l h_vmem=2G
#$ -l h_rt=720:00:00
#
#				2016-11-22
#				----------
#
# Here, the goal is to run ipyrad with the full dataset, and using a clustering
# threshold of 0.95. I will copy all required files to lluisvives:/scratch/jo/joiglu/
# and run it there, with qsub:
#
#   qsub -V -b n -pe openmpi 36 ./README.sh
#
# But, make sure you set the following variable to 'cluster' if running in the cluster:
WHEREAMIRUNNING=localhost
# WHEREAMIRUNNING=cluster

# Let's recall that the sequence data was produced with a GBS protocol, paired-end
# sequenced, and then merged. The GBS protocol included Y-shaped adapters. This
# contrasts with a less efficient GBS, where forward and reverse adapters are allowed
# to ligate randomly to either end of the fragments, and only fragments with both
# adapters get sequenced. In the latter case, it is clear that a fragment can be
# read from either of its complementary strands. Is this the case of the GBS with
# Y-shaped adapter?
#
# The fact is that a DNA fragment from GBS library prepared with Y-shaped adapters
# can also be read on either strand. To see this, I will represent below the genomic
# fragment strands as 'T' (top) and 'B' (bottom), to represent their complementarity.
# This double stranded fragment is ligated in both ends to double-stranded, Y-shaped
# adapters. Let's say that the top sequence of the adapter (5-prime to 3-prime) is
# '¿L', and the bottom is 'R!', where 'L' (for left) is complementary to R, but the
# question mark is not complementary to the exclamation mark. The double-stranded
# fragment with adapters ligated would look like this, with the strand on top always
# facing from 5-prime to 3-prime, and the one in the bottom, in reverse.
#
#				¿LTR!
#                               !RBL¿
#
# The amplification PCR uses the primers '{¡' and '[¿'. Upon denaturation, only primer
# '{¡' finds where to match during the first round.
#
#
#                                    ¿LTR!}     ___\     ¿LTR!}               1
#                                        ¡{        /     ?RBL¡{
#                                 /
#                                /
#        ¿LTR!    ___\    ¿LTR!}
#            ¡{      /    ?RBL¡{
#                                \
#                                 \
#                                  [¿           ___\    [¿LTR!}               2
#                                   ?RBL¡{         /    ]?RBL¡{
#
#
#
#
#                                  {¡LTR?       ___\    {¡LTR?]               3
#                                       ¿[         /    }!RBL¿[
#                                /
#                               /
#       {¡        ___\   {¡LTR?
#        !RBL¿       /   }!RBL¿
#                               \
#                                \
#                                  {¡           ___\    {¡LTR?                4
#                                  }!RBL¿          /    }!RBL¿
#
#
#
# Those are the four double-stranded fragments that result from PCR. The very top (1)
# and bottom (4) ones are shorter, since they miss one end or another. They stay in the
# reaction at low concentration because only one of their own strands keep producing them.
# The other one contributing to species (2) or (3). Species (2) and (3) are complete and
# reproduce by themselves along the PCR, becoming aboundant. Note that they are different.
# The sequencing primer for first read is something like '¿', binding to '?' and reading
# the top strand in species (2), but the bottom one in species 3.
#
# Thus, the data type should be GBS, which clusters the reads using reverse-complementary
# searches as well.
#

CURDIR=`pwd`
FASTQDIR=`pwd | sed 's/2016-11-22/2016-07-06b/'`
SAMPLE=(PipFe1 PipFe2 PipFe3 PipFe6 PipMa4 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5 Mol01 Mol02 Mol03 Mol04 Mol05)

if [ ! -d culex ]; then mkdir culex; fi
if [ ! -d culex/culex_fastqs ]; then mkdir culex/culex_fastqs; fi
for i in ${SAMPLE[@]}; do
   if [ ! -e culex/culex_fastqs/$i.fastq ]; then
      ln -s $FASTQDIR/$i'_setrimmed.fastq' culex/culex_fastqs/$i.fastq
   fi
done

if [ ! -e params-culex.txt ]; then
   ipyrad -n culex
   echo >> params-culex.txt
#                    --------ipyrad params file (v.0.4.1)--------------------------------------------
#                      culex                          ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
   sed -i "/## \[1\]/c culex                          ## [1] [project_dir]: Project dir (made in curdir if not present)"           params-culex.txt
                                                      ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                                                      ## [3] [barcodes_path]: Location of barcodes file
   sed -i "/## \[4\]/c culex/culex_fastqs/*.fastq     ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files"    params-culex.txt
#                      denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
                                                      ## [6] [reference_sequence]: Location of reference sequence file
   sed -i "/## \[7\]/c gbs                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc."               params-culex.txt
   sed -i "/## \[8\]/c CATG                           ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)" params-culex.txt
   sed -i "/## \[9\]/c 5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read"    params-culex.txt
   sed -i "/## \[10\]/c 30                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and...)" params-culex.txt
   sed -i "/## \[11\]/c 6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling"     params-culex.txt
   sed -i "/## \[12\]/c 1                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling"       params-culex.txt
   sed -i "/## \[13\]/c 100                            ## [13] [maxdepth]: Max cluster depth within samples"                       params-culex.txt
   sed -i "/## \[14\]/c 0.95                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly"       params-culex.txt
#                       0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
   sed -i "/## \[16\]/c 0                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)"   params-culex.txt
   sed -i "/## \[17\]/c 35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim"      params-culex.txt
   sed -i "/## \[18\]/c 2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences" params-culex.txt
   sed -i "/## \[19\]/c 5                              ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)"   params-culex.txt
   sed -i "/## \[20\]/c 10                             ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)"     params-culex.txt
   sed -i "/## \[21\]/c 1                              ## [21] [min_samples_locus]: Min # samples per locus for output"            params-culex.txt
   sed -i "/## \[22\]/c 40                             ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)"                    params-culex.txt
   sed -i "/## \[23\]/c 8                              ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)"             params-culex.txt
   sed -i "/## \[24\]/c 0.9                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)" params-culex.txt
   sed -i "/## \[25\]/c 0, 0                           ## [25] [edit_cutsites]: Edit cut-sites (R1, R2) (see docs)"                params-culex.txt
   sed -i "/## \[26\]/c 4, 4, 4, 4                     ## [26] [trim_overhang]: Trim overhang (see docs) (R1>, <R1, R2>, <R2)"     params-culex.txt
   sed -i "/## \[27\]/c *                              ## [27] [output_formats]: Output formats (see docs)"                        params-culex.txt
   sed -i "/## \[28\]/c $CURDIR/populations.txt        ## [28] [pop_assign_file]: Path to population assignment file"              params-culex.txt
fi

if [ ! -e populations.txt ]; then
   for i in `seq 0 16`; do
      echo ${SAMPLE[$i]} ${SAMPLE[$i]:0:3} >> populations.txt
   done
fi

touch checkpoints
if ! grep -q running checkpoints; then
   echo running >> checkpoints
   if [ $WHEREAMIRUNNING == 'cluster' ]; then
      ipyrad -p params-culex.txt -c 36 --MPI -s 234567
   fi
   if [ $WHEREAMIRUNNING == 'localhost' ]; then
#     echo "es una prova"
      ipyrad -p params-culex.txt -s 1234567
   fi
fi
