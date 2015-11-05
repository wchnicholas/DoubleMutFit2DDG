#R code
t <- read.table('data/PropInfo',header=1)
t$SS <- sub('E','purple',t$SS) #Beta Strand
t$SS <- sub('H','yellow',t$SS) #Alpha Helix
t$SS <- sub('T',colors()[240],t$SS) #Turn
t$SS <- sub('C',colors()[240],t$SS) #Coil
png('graph/PropSolSS.png')
barplot(rev(t$Sol),horiz=TRUE,col=rev(t$SS),names.arg=as.numeric(rev(t$Pos))+1,las=2,cex.names=0.6)
dev.off()
