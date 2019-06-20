source("/data/joiglu/bin/treemix-1.13/src/plotting_funcs.R")
# plot tree
plot_tree("culex_treemix_out/culex_stem")
dev.off()
file.rename(from = "Rplots.pdf", to = "culex_treemix_out/culex_tree.pdf")
# plot residues
plot_resid("culex_treemix_out/culex_stem", "poporder")
dev.off()
file.rename(from = "Rplots.pdf", to = "culex_treemix_out/culex_residues.pdf")
