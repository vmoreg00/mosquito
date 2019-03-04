#' Script that performs the ABBA-BABA test over the SNPs of Culex pipens
#' @author Moreno-González, V.
#' @date 2019-03-01

rm(list = ls())
source("/data/victor/mosquito/results/2019-02-23/src/ABBA_BABA_funs.R")

# Input files
wd <- "/data/victor/mosquito/results/2019-02-23/"
alt_freqs <- paste0(wd, "culex_alt_freqs.tsv")
chr_lengths <- paste0(wd, "culex_chr_lengths.tsv")

# Output files
abba_baba_out <- paste0(wd, "abba_baba.txt")

# Ploidies of the populations
pops <- c("molA1", "pipA4", "molM1", "torM2", "molS1", "mixS2", "mixS3", "torM4")
ploidies <- setNames(c(41, 41, 52, 56, 30, 26, 41, 41),
                     pops)

# Load the alternate allele frequencies table
freqs <- read.table(alt_freqs, sep = "\t", header = T, as.is = T)
names(freqs) <- c("chrom", "pos", pops)
cat("\nThere are", nrow(freqs), "SNPs in at least one population\n")

# join mix populations
freqs$mixS2.S3 <- apply(freqs[3:ncol(freqs)], 1, function(x){
  weighted.mean(x = c(x["mixS2"], x["mixS3"]),
                w = c(ploidies["mixS2"], ploidies["mixS3"]),
                na.rm = T)
})
freqs$mixS2.S3[is.nan(freqs$mixS2.S3)] <- NA
freqs$mixS2 <- NULL
freqs$mixS3 <- NULL

# Remove tor population (by definition, derive allele frequency in p0 is 0)
freqs$torM2 <- NULL
freqs$torM4 <- NULL

# filter (keep only positions where all populations have a value)
keep <- apply(freqs, 1, function(x) !any(is.na(x))) 
freqs <- freqs[keep, ]
cat(nrow(freqs), "sites remain after filtering")

# Chrom table:
chrom_table = read.table(chr_lengths)
chrom_lengths = chrom_table[,2]
names(chrom_lengths) = chrom_table[,1]

# Case 1: ((([mixS2+mixS3], molM1), pipA4), [torM4+torM2])
D1 <- abba.baba.test(freqs, "mixS2.S3", "molM1", "pipA4", 1e6, chrom_lengths)

# Case 2: ((([mixS2+mixS3], molM1), molA1), [torM4+torM2])
D2 <- abba.baba.test(freqs, "mixS2.S3", "molM1", "molA1", 1e6, chrom_lengths)

# Case 3: (((molS1, molM1), molA1), [torM4+torM2])
D3 <- abba.baba.test(freqs, "molS1", "molM1", "molA1", 1e6, chrom_lengths)

# Case 4: (((molS1, molM1), pipA4), [torM4+torM2])
D4 <- abba.baba.test(freqs, "molS1", "molM1", "pipA4", 1e6, chrom_lengths)

# Output
out_table <- rbind(D1, D2, D3, D4)
rownames(out_table) <- c("case1", "case2", "case3", "case4")
write.table(out_table, file = abba_baba_out, sep = "\t", quote = F)
