#!/bin/bash
#
#                               12-07-16
#                               ----------
#
#Para poder visualizar Snps y otras variaciones
#necesito convertir el formato BAM a VCF
LISTA=(Fe1r1 Fe2r1 Fe3r1 Fe4r1 Fe5r1 Fe6r1 Ma1r1 Ma2r1 Ma3r1 Ma4r1 Ma5r1 Ma6r1 Mol01r1 Mol02r1 Mol03r1 Mol04r1 Mol05r1 \
       Fe1r2 Fe2r2 Fe3r2 Fe4r2 Fe5r2 Fe6r2 Ma1r2 Ma2r2 Ma3r2 Ma4r2 Ma5r2 Ma6r2 Mol01r2 Mol02r2 Mol03r2 Mol04r2 Mol05r2)
#creo directorio
if [ ! -d VCF ];then mkdir VCF; fi
#Indexar la referencia para SAMtools
#primero me descargo la referencia
#despues si no existe ya la referencia indexada, lo indexo con samtools faidx
if [ ! -e reference.fa.fai ]; then
        wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/culex_quinquefasciatus/dna/Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz
	mv Culex_quinquefasciatus.CpipJ2.31.dna.toplevel.fa.gz reference.fa.gz
	samtools faidx reference.fa.gz 
fi
#conversion a vcf
for i in 0 1 2 3 4 5;do
	samtools mpileup -uf reference.fa  /home/student/Documents/mosquito/results/08-07-16/BAM/${LISTA[$i]}.bam > VCF/${LISTA[$i]}.mpileup
	#bcftools view -cg file.mpileup > file.vcf
done
