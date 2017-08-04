#!/bin/bash
#
#				2017-08-04
#				----------
#
# The problem of mapping genomic reads, with multiple origins, on a limited
# subset of loci is the low prior probability of a read mapping to the selected
# loci. This turns into a high false positive rate. In order to increase the
# probability of reads mapping to the selected loci, I should discard the reads
# mapping well elsewhere. I can blast the selected loci (2017-05-30/consensus.fasta)
# to the reference genome of C. quinquefasciatus in order to identify the positions
# of those loci. Then, I will map the reads, with global settings, to that reference
# genome, and select the ones that map around the loci identified, together with the
# ones that did not map anywhere, just in case.
