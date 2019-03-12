#' Functions used in the estimation of the standard error by the Jackknife
#' resampling method.
#' @autor Simon Martin
#' @modifyied Moreno-Gonzalez, V.
#' The modification introduced changes in the get_genome_blocks and 
#' get_jackknife_sd functions.
#' Now, all contigs which are smaller than the block_size will be ommited
#' to avoid an underestimation of the standard error.
get_genome_blocks <- function(block_size, chrom_lengths){
  block_starts <- sapply(chrom_lengths[chrom_lengths >= block_size], function(l){
    seq(1, l, block_size)
  })
  genome_blocks <- data.frame(start = unlist(block_starts),
                              end = unlist(block_starts) + block_size - 1,
                              chrom = rep(names(block_starts),
                                      sapply(block_starts, length)))
  genome_blocks <- genome_blocks[genome_blocks$end - genome_blocks$start == block_size - 1,]
  genome_blocks
}

get_genome_jackknife_indices <- function(chromosome, position, block_info){
  lapply(1:nrow(block_info), function(x){
    !(chromosome == block_info$chrom[x]
      & position >= block_info$start[x]
      & position <= block_info$end[x])
  })
}

get_jackknife_sd <- function(FUN, input_dataframe, jackknife_indices,
                             p1 = NA, p2 = NA, p3 = NA){
  n_blocks <- length(jackknife_indices)
  if(any(is.na(c(p1, p2, p3)))){
    overall_mean <- FUN(input_dataframe)
    blocks_means <- sapply(1:n_blocks, function(i){
      overall_mean * n_blocks - FUN(input_dataframe[jackknife_indices[[i]],]) * (n_blocks-1)
    })
  } else {
    overall_mean <- FUN(input_dataframe, p1, p2, p3)
    blocks_means <- sapply(1:n_blocks, function(i){
      overall_mean * n_blocks - FUN(input_dataframe[jackknife_indices[[i]],], p1, p2, p3) * (n_blocks-1)
    })
  }
  if(class(blocks_means) == "list"){
    sd(unlist(blocks_means))
  } else {
    sd(blocks_means)
  }
}
