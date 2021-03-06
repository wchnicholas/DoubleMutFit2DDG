#R code
t <- read.table('data/SingleECor',header=1)
png('graph/AllBGCor.png',res=100,height=300,width=300)
par(mar=c(2,2,0.3,0.3))
hist(t$Correlation,xlab='',ylab='',main='',col='grey')
dev.off()
png('graph/AllBGCorvsSASA.png',res=100,height=300,width=300)
par(mar=c(2,2,0.3,0.3))
plot(t$SASA,t$Correlation,cex=0.5)
dev.off()
png('graph/AllBGCorvsFit.png',res=100,height=300,width=300)
par(mar=c(2,2,0.3,0.3))
plot(t$Fitness,t$Correlation,cex=0.5)
dev.off()
