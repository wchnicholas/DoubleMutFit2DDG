#R code
png('/home/wchnicholas/Desktop/Manuscript/Stab/Fig/Fig/Model',res=100,height=300,width=500)
par(mar=c(0.2,0.2,0.3,0.3))
s <- seq(0,10,0.1)
plot(s,1/(1+exp(-(s-5)*1.8)),type='l',lwd=3,axes=F)
box(lwd=2)
dev.off()
