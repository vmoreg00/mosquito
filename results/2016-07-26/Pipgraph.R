clusters <- read.table('summary_clust.txt', header=TRUE)
attach(clusters)

png(filename='pipiens_clusters.png')

plot(Wclust[muestra == 'PipFe1_se'], total[muestra == 'PipFe1_se'], type='l')
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



