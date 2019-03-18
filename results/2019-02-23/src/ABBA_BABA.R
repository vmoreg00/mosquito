#' Script that performs the ABBA-BABA test over the SNPs of Culex pipens
#' @author Moreno-González, V.
#' @date 2019-03-01

rm(list = ls())
source("/data/victor/mosquito/results/2019-02-23/src/ABBA_BABA_funs.R")

############################## Initial varaibles #############################
# Input files
wd <- "/data/victor/mosquito/results/2019-02-23/"
alt_freqs <- paste0(wd, "culex_alt_freqs.tsv")
chr_lengths <- paste0(wd, "culex_chr_lengths.tsv")

# Output files
# abba_baba_out <- paste0(wd, "abba_baba.txt")

# Ploidies of the populations. Although during the execution of freebayes and
# CRISP (SNP callers), the ploidies of some population were reduced to
# simplify the computation requirements, here I indicate the original ploidies.
pops <- c("molA1", "pipA4", "molM1", "torM2", "molS1", "mixS2", "mixS3", "torM4")
ploidies <- setNames(c(224*2, 132*2, 26*2, 28*2, 15*2, 13*2, 64*2, 195*2),
                     pops)

######################### Allele frequencies table ###########################
# Load the alternate allele frequencies table
freqs <- read.table(alt_freqs, sep = "\t", header = T, as.is = T)
names(freqs) <- c("chrom", "pos", pops)
cat("\nThere are", nrow(freqs), "SNPs in at least one population\n")

# Join some populations to wiight their derived-allele frequencies
## "mix" populations
freqs$mixS2.S3 <- apply(freqs[3:ncol(freqs)], 1, function(x){
  weighted.mean(x = c(x["mixS2"], x["mixS3"]),
                w = c(ploidies["mixS2"], ploidies["mixS3"]),
                na.rm = T)
})
freqs$mixS2.S3[is.nan(freqs$mixS2.S3)] <- NA

## "tor" populations
freqs$torM2.M4 <- apply(freqs[3:ncol(freqs)], 1, function(x){
  weighted.mean(x = c(x["torM2"], x["torM4"]),
                w = c(ploidies["torM2"], ploidies["torM4"]),
                na.rm = T)
})
freqs$torM2.M4[is.nan(freqs$torM2.M4)] <- NA

# Filter: keep only positions where all used populations have a value != NA
keep <- apply(freqs, 1, function(x){
  !any(is.na(c(x["molA1"], x["pipA4"], x["molM1"], x["molS1"], x["mixS2.S3"])))
})
freqs_filtered <- freqs[keep, ]
cat(sum(keep == FALSE), "sites have been removed\n")
cat(nrow(freqs_filtered), "sites remain after filtering\n")

# Save
write.table(freqs_filtered,
            file = gsub(".tsv", "_filtered.tsv", alt_freqs),
            quote = F, sep = "\t")

########################### Chromosome lengths table #########################
# Chrom table:
chrom_table <- read.table(chr_lengths, col.names = c("chr", "len"))
chrom_lengths <- setNames(chrom_table$len, chrom_table$chr)

############################## ABBA_BABA Tests ###############################
# The same case is computed with different block sizes to have different
# jackknife estimates (more conservative while bigger blocks).
for(i in c(1e5, 5e5, 1e6)){
  ## Case 1: ((([mixS2+mixS3], molM1), pipA4), p0)
  cat("\n\nEstimating D and f under case 1 with block_size ==", i, ":\n")
  cat("\t((([mixS2+mixS3], molM1), pipA4), p0)\n")
  D1 <- abba.baba.test(freqs_filtered, "mixS2.S3", "molM1", "pipA4", i, chrom_lengths)
  print(D1)

  ## Case 2: ((([mixS2+mixS3], molM1), molA1), p0)
  cat("\n\nEstimating D and f under case 2 with block_size ==", i, ":\n")
  cat("\t((([mixS2+mixS3], molM1), molA1), p0)\n")
  D2 <- abba.baba.test(freqs_filtered, "mixS2.S3", "molM1", "molA1", i, chrom_lengths)
  print(D2)

  ## Case 3: (((molS1, molM1), pipA4), p0)
  cat("\n\nEstimating D and f under case 3 with block_size ==", i, ":\n")
  cat("\t(((molS1, molM1), pipA4), p0)\n")
  D3 <- abba.baba.test(freqs_filtered, "molS1", "molM1", "pipA4", i, chrom_lengths)
  print(D3)

  ## Case 4: (((molS1, molM1), molA1), p0)
  cat("\n\nEstimating D and f under case 4 with block_size ==", i, ":\n")
  cat("\t((([mixS2+mixS3], molM1), pipA4), (((molS1, molM1), molA1), p0)\n")
  D4 <- abba.baba.test(freqs_filtered, "molS1", "molM1", "molA1", i, chrom_lengths)
  print(D4)

  ## Output
  out_table <- rbind(D1, D2, D3, D4)
  rownames(out_table) <- c("case1", "case2", "case3", "case4")
  write.table(out_table, file = paste0(wd, "abba_baba_", i, ".txt"), sep = "\t", quote = F)
}
