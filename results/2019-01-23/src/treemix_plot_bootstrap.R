# Date: 2019-02-22
# Author: Moreno-Gonzalez, V.

library(ape)

# Variables
ml.tree <- "culex_treemix_out/culex_stem.treeout.gz"
bootstrap.dir <- "culex_treemix_bootstrap"
outdir <- "culex_treemix_bootstrap"

# Read the maximum likelihood tree
culex.tree <- read.tree(file = gzfile(ml.tree, "rt"))
culex.tree <- root(culex.tree, outgroup = c("4.torM2", "8.torM4"))

# Read the bootstrapped trees
boot.trees <- list()
for(i in list.files(bootstrap.dir, pattern = ".gz$", full.names = TRUE)){
  gz <- gzfile(i, "rt")
  boot.trees[[i]] <- root(read.tree(file = gz),
                     outgroup = c("4.torM2", "8.torM4"))
  close(gz)
}; rm(i)

# Plot the ML tree with the support on the edges
boot <- prop.clades(culex.tree, boot.trees)
pdf(file = paste0(outdir, "/culex_bootstrap_support.pdf"))
par(mar = rep(2, 4))
plot(culex.tree, main = "Culex pipiens treemix suport with bootstrap")
drawSupportOnEdges(boot)
dev.off()
