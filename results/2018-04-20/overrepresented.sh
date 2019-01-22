#!/bin/bash

# According to the fastQC report, in the file SRR2029628_1.fastq, there are
# some sequences overrepresented that do not correspond with any Illumina
# adapter:
# 	AAGAATTT
# 	AAAAATAC
# 	CTTAAATT
#
# Obviously, if I run trimmomatic, they are not going to dissapear.
#
# I have check if that sequences are also in the other files. They could be in
# a low frequency in the other files and FastQC cannot detect them. To do this,
# I have search the pattern of those sequences 4 or more times in the *_1.fastq
# files, and the reverse complementary in the *_2.fastq files.

DATA=/data/victor/mosquito/data/asgharian/raw

touch overrepseq.txt
for i in `seq 27 34`; do
      echo "$i\_forward" >> overrepseq.txt;
      egrep -c "(AAGAATTT){4,15}" $DATA/SRR20296$i\_1.fastq >> overrepseq.txt;
      egrep -c "(AAAAATAC){4,15}" $DATA/SRR20296$i\_1.fastq >> overrepseq.txt;
      egrep -c "(CTTAAATT){4,15}" $DATA/SRR20296$i\_1.fastq >> overrepseq.txt;
      echo "$i\_reverse" >> overrepseq.txt;
      egrep -c "(AAATTCTT){4,15}" $DATA/SRR20296$i\_2.fastq >> overrepseq.txt;
      egrep -c "(GTATTTTT){4,15}" $DATA/SRR20296$i\_2.fastq >> overrepseq.txt;
      egrep -c "(AATTTAAG){4,15}" $DATA/SRR20296$i\_2.fastq >> overrepseq.txt;
done

# After checking the file, I realized that that sequences (repeated 4 or more
# times) are in all the files!
# I have look for them in the NCBI, but there is no coincidence

# ------|---------------|---------------|---------
# 	| AAGAATTT	| AAAAATAC	| CTTAAATT
# ------|---------------|---------------|---------
# 27(+)	| 1336		| 4314		| 498
# 27(-)	| 1254		| 4153		| 644
# ------------------------------------------------
# 28(+)	| 242863	| 226928	| 244793
# 28(-)	| 233111	| 217456	| 249416
# ------------------------------------------------
# 29(+)	| 2928		| 10021		| 1076
# 29(-)	| 3243		| 7864		| 1187
# ------------------------------------------------
# 30(+)	| 26067		| 9971		| 1680
# 30(-)	| 26026		| 7251		| 3912
# ------------------------------------------------
# 31(+)	| 3101		| 24924		| 1155
# 31(-)	| 2859		| 22518		| 1452
# ------------------------------------------------
# 32(+)	| 1320		| 2440		| 342
# 32(-)	| 1530		| 1635		| 524
# ------------------------------------------------
# 33(+)	| 385		| 1840		| 274
# 33(-)	| 430		| 1775		| 322
# ------------------------------------------------
# 34(+)	| 3980		| 31996		| 1718
# 34(-)	| 3879		| 30159		| 2125
# ------------------------------------------------

# Probably, the sequences are artifacts or maybe, they are real. I will
# keep them unless they are problematic in the mapping or in the SNP
# detection.
