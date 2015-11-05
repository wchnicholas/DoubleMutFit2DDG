#R code
t <- read.table('data/PairwiseCorvsRSA',header=1)
png('graph/PairwiseCorvsRSA.png',res=100,height=600,width=600)
plot(t$RSA,t$Cor,pch=20,cex=0.05,xlab='RSA',ylab='Cor',main='')
dev.off()

