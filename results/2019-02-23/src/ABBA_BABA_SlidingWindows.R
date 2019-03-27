#' Script that performs the ABBA-BABA test over the SNPs of Culex pipens in
#' a sliding windows analysis
#' @author Moreno-González, V.
#' @date 2019-03-18

rm(list = ls())

library(ggplot2)
library(parallel)

source("src/ABBA_BABA_funs.R")

# Initial variables
input_tsv <- "culex_alt_freqs_filtered.tsv"
input_chroms <- "culex_chr_lengths.tsv"
window_size <- 1e+05
overlap <- 0
min_data <- 100
step <- window_size - overlap
cores <- 50

# Load the alternate allele frequencies table
freqs <- read.table(input_tsv, sep = "\t", header = T, as.is = T)

# Load chrom table
chrom_table <- read.table(input_chroms)
chrom_lengths <- chrom_table[,2]
names(chrom_lengths) <- chrom_table[,1]
chrom_table$cumm.length <- cumsum(chrom_table[, 2])

# Sliding window analysis (concurrent)
cl <- makeCluster(cores)
clusterExport(cl = cl, varlist = c("window_size", "overlap", "min_data", "step",
                                   "chrom_table", "chrom_lengths", "freqs", "f.stat",
                                   "abba.baba.test", "D.stat", "abba", "baba"))
D.sliding <- t(parallel::parSapply(
  cl = cl,
  X = as.character(chrom_table$V1),
  function(x){
    # If the contig is shorter than 5 * window_size, it is discarded
    if(any(chrom_lengths[x] < 5 * window_size, # Short contig
           sum(freqs$chrom == x) < 100)){      # Few or no SNPs in contig
      return(rep(NULL, 10))
    } else {
      windows_begin <- seq(1, chrom_lengths[x], by = window_size)
      D <- data.frame()
      k <- 1
      while(k < length(windows_begin)){
        # Start position of the window in freqs
        i <- windows_begin[k]
        freq_begin <- which(freqs$chrom == x & freqs$pos >= i)[1]
        # Last position of the window in freqs
        j <- i + window_size
        freq_end <- tail(which(freqs$chrom == x & freqs$pos <= j), 1)
        # Number of SNPs
        n.snp <- freq_end - freq_begin
        if(n.snp < 100){
          # Extention of the window and re-computation of the windows' begin
          freq_end <- freq_begin + 100
          offset <- freqs$pos[freq_end] - j
          windows_begin <- windows_begin + offset
        }
        if(j > chrom_lengths[x]){
          D <- rbind(D, rep(NA, 10))
        } else {
          # Execute ABBA_BABA test in 2 cases
          D_case1 <- c(x, i, abba.baba.test(freqs[freq_begin:freq_end,],
                                            "mixS2.S3", "molM1", "pipA4",
                                            sliding = TRUE))
          D_case2 <- abba.baba.test(freqs[freq_begin:freq_end,],
                                    "mixS2.S3", "molM1", "molA1",
                                    sliding = TRUE)
          D <- rbind(D, c(D_case1, D_case2), stringsAsFactors = F)
          colnames(D) <- c("chr", "pos",
                           "D_1", "f_1", "ABBA_1", "BABA_1",
                           "D_2", "f_2", "ABBA_2", "BABA_2")
        }
        k <- k + 1
      }
      return(D)
    }
  }))
stopCluster(cl)

# Format the table and save it
D.sliding1 <- data.frame()
for(i in D.sliding){
  D.sliding1 <- rbind(D.sliding1, i)
}
for(i in 2:ncol(D.sliding1)){
  D.sliding1[, i] <- as.numeric(D.sliding1[, i])
}
## Differences between the f statistic in both cases
D.sliding1$diff_f <- D.sliding1$f_1 - D.sliding1$f_2
## "Genomic position" if contigs were sorted. For representation purposes.
## The offset of the extended windows is very low, so it can be asumed that
## all the windows are 100000bp length.
D.sliding1$pos_genomic <- seq(from = 1, by = 100000,
                              length.out = nrow(D.sliding1))
## Save
write.table(D.sliding1, file = "abba_baba_slidingWindows.tsv", sep = "\t",
            quote = F)

# Plot the results
facets <- 5
D.sliding1$facet <- cut(D.sliding1$pos_genomic, facets, seq(1, facets))

g <- ggplot(data = D.sliding1, aes(x = pos_genomic)) +
  geom_line(aes(y = diff_f)) +
  facet_wrap(vars(facet), ncol = 1, scales = "free_x",
             strip.position = "right") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE),
                     name = "Genome position") +
  scale_y_continuous(name = bquote(~f[1] - ~f[2]), limits = c(-1, 1)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank())
ggsave(filename = "abba_baba_slidingWindows.png", plot = g, device = "png",
       width = 250, height = 175, units = "mm", dpi = 333)
