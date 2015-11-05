#R code
#Plot the histogram for the R_literature for the sim result
histplot <- function(t,break1,break2){
  t <- t[which(sol$k >= 10),]
  hist(t$X0.to.1,breaks=break1,xlim=c(-0.2,1),main='',xlab='')
  abline(v=0.52,lty=2)
  hist(t$X0.4.to.0.8,breaks=break2,xlim=c(-0.2,1),main='',xlab='')
  abline(v=0.52,lty=2)
  }

png('graph/CorSimResult_nolabel.png',res=100,height=250,width=400)
par(mfrow=c(4,1),mar=c(2.2,4.2,0.2,0.2))
sol <- read.table('data/SimResult_Sol',sep="\t",header=1)
histplot(sol,20,20)
hyd <- read.table('data/SimResult_Hydro',sep="\t",header=1)
histplot(hyd,200,50)
dev.off()
