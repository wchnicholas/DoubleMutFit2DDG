#R code

plothist <- function(infile,outfile,breaks){
  t <- read.table(infile,header=1)
  png(outfile,res=100,height=100,width=250)
  par(mar=c(2,2,0.3,0.3))
  hist(t[,3],xlab='',ylab='',main='',col='grey',xlim=c(0,1),breaks=breaks)
  dev.off()
  }

plothist('data/GroupIECor','graph/GroupISASA.png',20)
plothist('data/GroupIIECor','graph/GroupIISASA.png',20)
plothist('data/GroupIIIECor','graph/GroupIIISASA.png',20)

g1 <- read.table('data/GroupIECor')
g2 <- read.table('data/GroupIIECor')
g3 <- read.table('data/GroupIIIECor')
