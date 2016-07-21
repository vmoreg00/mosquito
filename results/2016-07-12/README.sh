#!/bin/bash
#
#                               12-07-16
#                               ----------
#
#PATHS
DIR=../2016-07-08
#Para poder visualizar Snps y otras variaciones
#necesito convertir el formato BAM a VCF
LISTA=(PipFe1 PipFe2 PipFe3 PipFe4 PipFe5 PipFe6 PipMa1 PipMa2 PipMa3 PipMa4 PipMa5 PipMa6 Mol01 Mol02 Mol03 Mol04 Mol05)

#Indexar la referencia para SAMtools
#primero me descargo la referencia
#despues si no existe ya la referencia indexada, lo indexo con samtools faidx
if [ ! -e reference.fa.fai ]; then
        wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz
	mv Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz reference.fa.gz
	gunzip reference.fa.gz
	samtools faidx reference.fa
fi
#Debo tener indexadas y ordenadas las muestras por posicion y cromosomas
for i in 0 1 2 3 4 5; do
	if [ ! -e  ${LISTA[$i]}'_sorted' ]; then
		samtools sort  $DIR/${LISTA[$i]}.bam  ${LISTA[$i]}'_sorted'
	fi
	if [ ! -e ${LISTA[$i]}'_sorted'.bam.bai ]; then
		samtools index   ${LISTA[$i]}'_sorted'.bam
	fi
done
#conversion a vcf
if [ ! -e mosquito.vcf ]; then
	samtools mpileup -gf reference.fa  *_sorted.bam > mosquito.mpileup
	bcftools view -Acg mosquito.mpileup > mosquito.vcf
fi
