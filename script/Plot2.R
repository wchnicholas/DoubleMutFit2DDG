#R code
t <- read.table('data/SingleECor',header=1)
png('graph/AllBGsSubset.png',res=100,height=150,width=600)
par(mar=c(2,2,0.3,0.3))
stripchart(t[which(t$SASA < 0.03 & t$Fitness > 0.4 & t$Fitness < 0.6),]$Correlation,xlab='',ylab='',main='',pch=20,method='jitter',jitter=0.01)
dev.off()
