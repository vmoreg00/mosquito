#' Set of functions used in the ABBA BABA test. This code is
#' inspired in the workshop of evomics.org (Tutorial: ABBA-BABA statistics)
#' developed by Simon Martin
#' @author Moreno-Gonzalez, V.
#' @date 2019-03-01

source("/data/victor/mosquito/results/2019-02-23/src/jackknife.R")
#' @description Computation of the ABBA proportion
abba <- function(p1, p2, p3){
  (1 - p1) * p2 * p3
}

#' @description Computation of the BABA proportion
baba <- function(p1, p2, p3){
  p1 * (1 - p2) * p3
}

#' @description Computation of the D statistic based of the ABBA and BABA 
#' patterns
D.stat <- function(df){
  (sum(df$ABBA) - sum(df$BABA)) / (sum(df$ABBA) + sum(df$BABA))
}

#' @description Carries out the ABBA BABA test, computes de D statistics,
#' its deviation and signification based on the Jacknife procedure
abba.baba.test <- function(freqs, p1, p2, p3, block.size = 1e6, chr.lengths){
  # D statistic
  ABBA_BABA <- data.frame(ABBA = abba(freqs[, p1], freqs[, p2], freqs[, p3]),
                          BABA = baba(freqs[, p1], freqs[, p2], freqs[, p3]))
  D <- D.stat(ABBA_BABA)
  
  # Jackknife
  blocks <- get_genome_blocks(block_size = block.size, 
                              chrom_lengths = chr.lengths)
  n_blocks <- nrow(blocks)
  indices <- get_genome_jackknife_indices(chromosome = freqs$chrom,
                                          position = freqs$pos,
                                          block_info = blocks)
  D_sd <- get_jackknife_sd(FUN = D.stat, 
                           input_dataframe = ABBA_BABA,
                           jackknife_indices = indices)
  D_err <- D_sd/sqrt(n_blocks)
  D_Z <- D / D_err
  D_p <- 2*pnorm(-abs(D_Z))
  
  # Output
  D_total <- setNames(c(D, D_sd, D_err, D_Z, D_p),
                      c("D", "D_sd", "D_err", "D_Z", "D_p"))
  return(D_total)
}

