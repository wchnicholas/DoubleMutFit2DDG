#R code
library(marray)
library(gplots)
t <- read.table('data/LandCorMatrix',header=1,row.names=1)
t <- t(t)
mcol <- maPalette(low="blue", mid="white", high="red")
png('graph/LandCorHeatmap.png',res=1200,height=7200,width=7200)
#png('graph/LandCorHeatmap.png',res=100,height=600,width=600)
#heatmap.2(t,col=mcol,dendrogram='both',trace="none",labCol=F,labRow=F,keysize=1,density.info='none')
#hm <- heatmap.2(t,col=mcol,dendrogram='both',trace="none",labCol=F,labRow=F,keysize=1,density.info='none',cexRow=0.1)
hm <- heatmap.2(t,col=mcol,dendrogram='both',trace="none",keysize=1,density.info='none',cexRow=0.05)
dev.off()
clus <- t[rev(hm$rowInd), hm$colInd]
write.table(clus,'data/LandhcMatrix')


#t <- read.table('data/MutPotCorMatrix',header=1,row.names=1)
#t <- t(t)
#mcol <- maPalette(low="blue", mid="white", high="red")
#png('graph/MutPotCorHeatmap.png',res=1200,height=7200,width=7200)
#hm <- heatmap.2(t,col=mcol,dendrogram='both',trace="none",labCol=F,keysize=1,density.info='none',cexRow=0.05)
#dev.off()
#clus <- t[rev(hm$rowInd), hm$colInd]
#write.table(clus,'data/MutPothcMatrix')
