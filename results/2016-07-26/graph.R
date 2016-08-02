clusters <- read.table('summary_clust.txt', header=TRUE)
attach(clusters)

png(filename='molestus_clusters.png')

plot(Wclust[muestra == 'Mol01_se'], total[muestra == 'Mol01_se'], type='l')
lines(Wclust[muestra == 'Mol02_se'], total[muestra == 'Mol02_se'])
lines(Wclust[muestra == 'Mol03_se'], total[muestra == 'Mol03_se'])
lines(Wclust[muestra == 'Mol04_se'], total[muestra == 'Mol04_se'])
lines(Wclust[muestra == 'Mol05_se'], total[muestra == 'Mol05_se'])

dev.off()
