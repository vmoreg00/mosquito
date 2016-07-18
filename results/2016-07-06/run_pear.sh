#!/bin/bash

pear  -f $1'_R1.fastq' \
      -r $1'_R2.fastq' \
      -o merged/$1 \
      -v 10 \
      -q 10 \
      -j 1 \
      --memory 2G

vsearch --fastx_revcomp merged/$1.unassembled.reverse.fastq \
        --fastqout merged/$1.unassembled.reverse.reversed.fastq

mv merged/$1.unassembled.reverse.reversed.fastq merged/$1.unassembled.reverse.fastq
