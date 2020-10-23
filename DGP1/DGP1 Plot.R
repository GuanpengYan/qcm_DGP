library(readxl)
library(tidyverse)
library(devEMF)

File="DGP1"
M=read.csv(file=paste(File,"Cover.csv"),row.names=1)
rowM=dim(M)[1]
rowP=3
colP=rowM/rowP
qrf.plot=function(){
  par(mfrow=c(rowP,colP))
  for (i in 1:rowM){
    Main=paste(File,", N=") %>% paste(.,i*10,sep="")
    plot(10*(1:9),M[i,],main =Main,xlab = "Pretreatment Period (T0)",ylab = "Coverage Probability",ylim=c(0.5,1),xaxt="n")
    axis(1, at = seq(10, 90, by = 10))
    lines(10*(1:9),M[i,])
    abline(h=0.95,lty=2)
    abline(h=0.9,lty=3)
  }
}
qrf.plot()
emf(file=paste(File,"Photo.emf"), width=6, height=8,bg="transparent")
qrf.plot()
dev.off()
