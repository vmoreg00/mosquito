# DESCRIPTION #################################################################
#
# Script to plot a dendrogram (as phylogeny) with the averaged Fst values
# obtained in 10Kb non-overlapping sliding windows with PoPoolation2.
#
# Author: Moreno-Gonzalez, V.
# Date: 2019-06-11

rm(list = ls())
library(ape)

# MAIN ########################################################################
## Load the data
labs <- c("C. pipiens f. molestus\n(molA1)", "C. pipiens f. pipiens\n(pipA4)",
          "C. pipiens f. molestus\n(molM1)", "C. torrentium\n(torM2)",
          "C. pipiens f. molestus\n(molS1)", "C. pipiens mixed\n(mixS2)",
          "C. pipiens mixed\n(mixS3)", "C. torrentium\n(torM4)")
fst <- read.table("culex_fst_mean.csv", header = 1, row.names = 1, sep = "\t")
colnames(fst) <- labs
rownames(fst) <- labs

## Make the cluster analysis (stored as a phylo object)
fst <- as.dist(t(fst))
fst_clust <- as.phylo(hclust(fst))
## Rotate a interior node for representation purposes
fst_clust <- rotate(phy = fst_clust, node = 12)

## Plot the phylogeny
png(filename = "culex_fst_phylo.png", 
    width = 20, height = 20, units = "cm", res = 333)
plot(fst_clust, use.edge.length = T, 
     edge.width = 2, cex = 1,
     main = paste0("Dendrogram according to ", 
                   "average Fst values in a 10Kb\n", 
                   "non-overlapping sliding windows"),
     cex.main = 0.8)
dev.off()


