lusters <- read.table('summary_clust.txt', header=TRUE)
attach(clusters)

png(filename='molestus_clusters.png')

plot(Wclust[muestra == 'Mol01_se'], total[muestra == 'Mol01_se'], type='l',
   xlab = "Clustering threshold (%)", ylab = "Number of clusters", main = 'Molestus samples')
lines(Wclust[muestra == 'Mol02_se'], total[muestra == 'Mol02_se'])
lines(Wclust[muestra == 'Mol03_se'], total[muestra == 'Mol03_se'])
lines(Wclust[muestra == 'Mol04_se'], total[muestra == 'Mol04_se'])
lines(Wclust[muestra == 'Mol05_se'], total[muestra == 'Mol05_se'])

dev.off()

png(filename='pipiens_clusters.png')

plot(Wclust[muestra == 'PipFe1_se'], total[muestra == 'PipFe1_se'], type='l',
   xlab = 'Clustering threshold (%)', ylab = 'Number of clusters', main = 'Pipiens samples')
lines(Wclust[muestra == 'PipFe2_se'], total[muestra == 'PipFe2_se'])
lines(Wclust[muestra == 'PipFe3_se'], total[muestra == 'PipFe3_se'])
lines(Wclust[muestra == 'PipFe4_se'], total[muestra == 'PipFe4_se'])
lines(Wclust[muestra == 'PipFe5_se'], total[muestra == 'PipFe5_se'])
lines(Wclust[muestra == 'PipFe6_se'], total[muestra == 'PipFe6_se'])
lines(Wclust[muestra == 'PipMa1_se'], total[muestra == 'PipMa1_se'])
lines(Wclust[muestra == 'PipMa2_se'], total[muestra == 'PipMa2_se'])
lines(Wclust[muestra == 'PipMa3_se'], total[muestra == 'PipMa3_se'])
lines(Wclust[muestra == 'PipMa4_se'], total[muestra == 'PipMa4_se'])
lines(Wclust[muestra == 'PipMa5_se'], total[muestra == 'PipMa5_se'])
lines(Wclust[muestra == 'PipMa6_se'], total[muestra == 'PipMa6_se'])


dev.off()
