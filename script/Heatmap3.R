#R code
library(marray)
library(gplots)
library(stringr)
t <- read.table('data/PropCorMatrix',header=1,row.names=1)
t <- t(t)
row.names(t) <- as.numeric(str_sub(row.names(t),2,-1))+1
mcol <- maPalette(low="blue", mid="white", high="red")
png('graph/PropCor.png',res=200,height=1200,width=1200)
hm <- heatmap.2(t,col=mcol,dendrogram='both',trace="none",labCol=F,keysize=1,density.info='none',cexRow=0.6)
dev.off()
clus <- t[rev(hm$rowInd), hm$colInd]
write.table(clus,'data/PropCorhc')
