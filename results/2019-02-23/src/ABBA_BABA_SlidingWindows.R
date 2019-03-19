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
                                   "chrom_table", "freqs", "f.stat",
                                   "abba.baba.test", "D.stat", "abba", "baba"))
D.sliding <- t(parallel::parSapply(
  cl = cl,
  X = seq(from = 1, by = step,
          length.out = round(max(chrom_table$cumm.length) / step, 0)),
  function(x){
    # Determine the start and end position of the window IN THE GENOME
    ## Start and end position of the window in the genome
    begin <- x
    end <- x + window_size
    ## First chromosome in window
    chr_begin <- which(chrom_table$cumm.length > begin)
    chr_begin <- as.character(chrom_table[, 1][chr_begin[1]])
    ## Start position in the genome (3rd column of chrom_table)
    if(begin == 1){
      pos_begin <- 1
    } else {
      pos_begin <- which(chrom_table[, 1] == chr_begin)
      pos_begin <- abs(chrom_table[pos_begin - 1, 3] - begin)
    }
    ### If chr_begin has no SNP (not in freqs), the next will be selected
    while(!(chr_begin %in% freqs$chrom)){
      chr_begin <- 
        as.character(chrom_table[which(chrom_table[, 1] == chr_begin) + 1, 1])
      pos_begin <- 1
    }
    ## Last chromosome in window
    chr_end <- which(chrom_table$cumm.length > end)
    chr_end <- as.character(chrom_table[, 1][chr_end[1]])
    ## End position in the genome (3rd column of chrom_table) 
    pos_end <- which(chrom_table[, 1] == chr_end)
    pos_end <- abs(chrom_table[pos_end - 1, 3] - end)
    ### If chr_end has no SNP (not in freqs), the previous will be selected
    while(!(chr_end %in% freqs$chrom)){
      chr_end <-
        as.character(chrom_table[which(chrom_table[, 1] == chr_end) - 1, 1])
      pos_end <- chrom_table[which(chrom_table[, 1] == chr_end), 2]
    }

    # Determine which SNPs are between start and end position in freqs_table
    ## First SNP in window
    freq_begin <- which(freqs$chrom == chr_begin & freqs$pos >= pos_begin)[1]
    ### If freq_begin does not exist, the first of the next chr is selected
    if(is.na(freq_begin)){
      freq_begin <- tail(which(freqs$chrom == chr_begin), 1) + 1
    }
    ## Last SNP in window
    freq_end <- which(freqs$chrom == chr_end & freqs$pos <= pos_end)[1]
    ### If freq_end does not exist, the last SNP of the previous chr is selected
    if(is.na(freq_end)){
      freq_end <- which(freqs$chrom == chr_end)[1] - 1
    }

    # ABBA BABA test in 2 cases
    if(freq_end - freq_begin > min_data){
      D_case1 <- c(x, abba.baba.test(freqs[freq_begin:freq_end,],
                                     "mixS2.S3", "molM1", "pipA4",
                                     sliding = TRUE))
      D_case2 <- abba.baba.test(freqs[freq_begin:freq_end,],
                                "mixS2.S3", "molM1", "molA1",
                                sliding = TRUE)
    } else {
      D_case1 <- c(x, 0, 0, 0, 0)
      D_case2 <- c(0, 0, 0, 0) # pos is ommited as it is the same as in case 1
    }
    D <- c(D_case1, D_case2)
    D
  }))
stopCluster(cl)

# Format the table and save it
D.sliding <- as.data.frame(D.sliding)
names(D.sliding) <- c("pos", "D_1", "f_1", "ABBA_1", "BABA_1", 
                      "D_2", "f_2", "ABBA_2", "BABA_2")
write.table(D.sliding, file = "abba_baba_slidingWindows.tsv", sep = "\t", 
            quote = F)

# Plot the results
facets <- 5
D.sliding$facet <- cut(D.sliding$pos, facets, seq(1, facets))

g <- ggplot(data = D.sliding, aes(x = pos)) +
  geom_line(aes(y = D_1, colour = "p3 = pipiens from Aleksin")) +
  geom_line(aes(y = D_2, colour = "p3 = molestus from Aleksin")) +
  facet_wrap(vars(facet), ncol = 1, scales = "free_x",
             strip.position = "right") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE),
                     name = "Genome position") +
  scale_y_continuous(name = "D statistic", limits = c(-1, 1)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank())
ggsave(filename = "abba_baba_slidingWindows.png", plot = g, device = "png",
       width = 250, height = 175, units = "mm", dpi = 333)
