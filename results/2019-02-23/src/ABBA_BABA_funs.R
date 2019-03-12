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

#' @description Computes the f statistic, which is an estimate of the
#' admixture proportion. Chosen implementation is the conservative form,
#' f_{hom} (see Martin et al., 2014; doi:10.1093/molbev/msu269) where it is
#' assumed that under complete introgression, the allele frequencies of
#' P2 and P3 will homogenize.
f.stat <- function(freqs, p1, p2, p3){
  # Observed difference (obs.S) between ABBA and BABA
  obs.abba <- abba(freqs[, p1], freqs[, p2], freqs[, p3])
  obs.baba <- baba(freqs[, p1], freqs[, p2], freqs[, p3])
  obs.S <- sum(obs.abba) - sum(obs.baba)

  # Expected difference (exp.S) between ABBA and BABA
  exp.abba <- abba(freqs[, p1], freqs[, p3], freqs[, p3])
  exp.baba <- baba(freqs[, p1], freqs[, p3], freqs[, p3])
  exp.S <- sum(exp.abba) - sum(exp.baba)

  # Estimation of f statistic
  f <- obs.S / exp.S
  return(f)
}

#' @description Carries out the ABBA BABA test, computes de D statistics,
#' its deviation and signification based on the Jacknife procedure
abba.baba.test <- function(freqs, p1, p2, p3, block.size = 1e6, chr.lengths){
  # D statistic
  ABBA_BABA <- data.frame(ABBA = abba(freqs[, p1], freqs[, p2], freqs[, p3]),
                          BABA = baba(freqs[, p1], freqs[, p2], freqs[, p3]))
  D <- D.stat(ABBA_BABA)

  # f statistic
  f <- f.stat(freqs, p1, p2, p3)
  
  # Jackknife...
  blocks <- get_genome_blocks(block_size = block.size, 
                              chrom_lengths = chr.lengths)
  n_blocks <- nrow(blocks)
  indices <- get_genome_jackknife_indices(chromosome = freqs$chrom,
                                          position = freqs$pos,
                                          block_info = blocks)
  ## ... for D statistic
  D_sd <- get_jackknife_sd(FUN = D.stat, 
                           input_dataframe = ABBA_BABA,
                           jackknife_indices = indices)
  D_err <- D_sd / sqrt(n_blocks)
  D_Z <- D / D_err
  D_p <- 2 * pnorm(-abs(D_Z))
  ## ... for f statistic
  f_sd <- get_jackknife_sd(FUN = f.stat,
                           p1 = p1, p2 = p2, p3 = p3,
                           input_dataframe = freqs,
                           jackknife_indices = indices)
  f_err <- f_sd / sqrt(n_blocks)
  f_Z <- f / f_err
  f_p <- 2 * pnorm(-abs(f_Z))

  # Output
  summary_total <- setNames(c(f, f_sd, f_err, f_Z, f_p,
                              D, D_sd, D_err, D_Z, D_p),
                            c("f", "f_sd", "f_err", "f_Z", "f_p",
                              "D", "D_sd", "D_err", "D_Z", "D_p"))
  return(summary_total)
}
