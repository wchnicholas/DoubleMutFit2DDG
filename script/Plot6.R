#R code
coloring <- function(cor){
  if (cor > 0.85){return (rgb(1,0,0,0.7))}
  else if (cor > 0.75){return (rgb(0,0,1,0.7))}
  else {return (rgb(0.5,0.5,0.5,0.2))}
  }

plotting <- function(filename,main){
  t <- read.table(filename,header=1)
  col <- mapply(coloring,t$Correlation)
  png(main,res=100,height=300,width=300)
  par(mar=c(4.2,4.2,0.2,0.2))
  plot(t$Fitness,t$SASA,col=col,cex=0.8,pch=20,xlim=c(0,1),ylim=c(0,1),xlab='Fitness',ylab='RSA')
  dev.off()
  }

plotting('data/SingleECor','graph/RSAvsFit_All.png')
plotting('data/GroupIECor','graph/RSAvsFit_GroupI.png')
plotting('data/GroupIIECor','graph/RSAvsFit_GroupII.png')
plotting('data/GroupIIIECor','graph/RSAvsFit_GroupIII.png')
