#!/bin/bash
#
#				2016-07-21
#				----------
# We will run pyrad to the assembled samples we got in the directory 2016-07-06b.
# Most of the reads had been assembled so we decided  to work with those files.
# The input file in pyrad is params.txt, so if we have not  it, we must create it

if [ ! -e params.txt ]; then
	pyrad -n
	sed -i "/## 2. /c\../2016-07-06b/*_trimmed.fastq 	## 2. path to raw data files" params.txt
	sed -i '/## 5. /c\muscle	    	## 5. path to call muscle' params.txt
	sed -i '/## 7. /c\6                 	## 7. N processors (parallel) (all)' params.txt
	sed -i '/## 8. /c\2                 	## 8. Mindepth: min coverage for a cluster (s4,s5)' params.txt
   	sed -i '/## 9. /c\4                 	## 9. maxN: max number of Ns in reads (s2)' params.txt
  	sed -i '/## 10. /c\.90              	## 10. Wclust: clustering threshold as a decimal (s3,s6)' params.txt
   	sed -i '/## 11. /c\rad              	## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merged (all)' params.txt
   	sed -i '/## 13. /c\5                	## 13. MaxSH: max inds with shared hetero site (s7)' params.txt
   	sed -i '/## 14. /c\mosq             	## 14. Prefix name for final output (no spaces) (s7)' params.txt
	sed -i "/## 18./c\../2016-07-06b/*_trimmed.fastq  ##18.opt.: loc. of de-multiplexed data (s2)" params.txt
	sed -i '/## 21. /c\1                	## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)' params.txt
   	sed -i '/## 30. /c\*                	## 30.opt.: Output formats... (s7)' params.txt
   	sed -i '/## 36. /c\1                	## 36.opt.: vsearch max. threads per job (def.=6; see docs) (s3,s6)' params.txt
	cd ..
fi
#Now we  execute pyrad
pyrad -p params.txt -s 234567 1> pyrad.log 2> pyrad.err

