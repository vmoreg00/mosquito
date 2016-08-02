#!/bin/bash
#
#                               2016-07-26
#                               ----------
# We will run pyrad to the assembled samples we got in the directory 2016-07-06b.
# Most of the reads had been assembled so we decided  to work with those files.
# The input file in pyrad is params.txt, so if we have not  it, we must create it

# Random sample of input

DIR=../2016-07-06b
LISTA=(PipFe1 PipFe2 PipFe3 PipFe6 PipMa4 PipFe4 PipMa3 PipMa1 PipMa2 PipMa5 PipMa6 PipFe5 Mol01 Mol02 Mol03 Mol04 Mol05)
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
	if [ ! -e  ${LISTA[$i]}'_se.fastq' ]; then
		sample_seqs -t fastq -o  ${LISTA[$i]}'_se.fastq' $DIR/${LISTA[$i]}'_setrimmed.fastq' -n 100000  
	fi
done
# Creation of params.txt
if [ ! -e params.txt ]; then
        pyrad -n
        sed -i "/## 2. /c\*_se.fastq      ## 2. path to raw data files" params.txt
        sed -i '/## 5. /c\muscle                ## 5. path to call muscle' params.txt   
        sed -i '/## 7. /c\6                     ## 7. N processors (parallel) (all)' params.txt
        sed -i '/## 8. /c\2                     ## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
        sed -i '/## 9. /c\4                     ## 9. maxN: max number of Ns in reads (s2)' params.txt
        sed -i '/## 10. /c\.70                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
        sed -i '/## 11. /c\rad                  ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
        sed -i '/## 13. /c\5                    ## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
        sed -i '/## 14. /c\mosq                 ## 14. Prefix name for final output (no spaces) (s7)' params.txt
        sed -i "/## 18./c\*_se.fastq  ##18.opt.: loc. of de-multiplexed data (s2)" params.txt
        sed -i '/## 21. /c\1                    ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
        sed -i '/## 30. /c\*                    ## 30.opt.: Output formats... (s7)' params.txt
        sed -i '/## 36. /c\1                    ## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
fi

# Calling step2
pyrad  -p params.txt -s2 1> pyrad.log 2> pyrad.err 
# Pyrad's execution each time with a different Wclust
# Wclust .70
sed -i '/## 10. /c\.70                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
if [ ! -e stats/s3.70.txt ]; then
	pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err 
	mv stats/s3.clusters.txt stats/s3.70.txt
fi
# Wclust .75
sed -i '/## 10. /c\.75                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
if [ ! -e stats/s3.75.txt ]; then
        pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err  
        mv stats/s3.clusters.txt stats/s3.75.txt
fi
# Wclust .80
sed -i '/## 10. /c\.80                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt 
if [ ! -e stats/s3.80.txt ]; then
        pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err  
	 mv stats/s3.clusters.txt stats/s3.80.txt
fi 
# Wclust .85
sed -i '/## 10. /c\.85                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
if [ ! -e stats/s3.85.txt ]; then
        pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err  
        mv stats/s3.clusters.txt stats/s3.85.txt
fi
# Wclust .90
sed -i '/## 10. /c\.90                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt 
if [ ! -e stats/s3.90.txt ]; then
        pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err  
        mv stats/s3.clusters.txt stats/s3.90.txt
fi
# Wclust .95
sed -i '/## 10. /c\.95                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt 
if [ ! -e stats/s3.95.txt ]; then
	pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err  
        mv stats/s3.clusters.txt stats/s3.95.txt
fi 
# Wclust .99
sed -i '/## 10. /c\.99                  ## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt 
if [ ! -e stats/s3.99.txt ]; then
        pyrad -p params.txt -s 3 1>> pyrad.log 2>> pyrad.err  
        mv stats/s3.clusters.txt stats/s3.99.txt
fi 

 After everything is done, we run the summary:
if [ ! -e summary_clust.txt ]; then
   if [ ! -e archivos.txt ]; then
	ls -1 stats/s3* > archivos.txt
  fi
   ./summary.py archivos.txt
   rm archivos.txt
fi


