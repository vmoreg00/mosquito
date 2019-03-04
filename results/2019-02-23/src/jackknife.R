get_genome_blocks <- function(block_size, chrom_lengths) {
    block_starts <- sapply(chrom_lengths, function(l) seq(1, l, block_size))
    data.frame(start = unlist(block_starts),
               end = unlist(block_starts) + block_size - 1,
               chrom = rep(names(block_starts), sapply(block_starts, length)))
    }

get_genome_jackknife_indices <- function(chromosome, position, block_info){
    lapply(1:nrow(block_info), function(x) !(chromosome == block_info$chrom[x] &
                                             position >= block_info$start[x] &
                                             position <= block_info$end[x]))
    }


get_jackknife_sd <- function(FUN, input_dataframe, jackknife_indices){
    n_blocks <- length(jackknife_indices)
    overall_mean <- FUN(input_dataframe)
    sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - FUN(input_dataframe[jackknife_indices[[i]],])*(n_blocks-1)))
    }
