rm(list = ls())
library(latex2exp)

## define plot data
xlim <- c(0,1)
ylim <- c(0,1)

px <- c(0.2,0.3,0.4)
py <- c(0,0,0)

lx.buf <- 5
lx<-c(0.1,0.3,0.5)
#lx <- seq(xlim[1]+lx.buf,xlim[2]-lx.buf,len=length(px))
ly <- 0.2

## create basic plot outline
par(xaxs='i',yaxs='i',mar=c(5,1,1,1))
plot(NA,xlim=xlim,ylim=ylim,axes=F,ann=F)
axis(1)

## plot elements
segments(px,py,lx,ly)
points(px,py,pch=16,xpd=NA)
text(lx,ly,labels=c(TeX('$\\theta_0$'),TeX('$\\theta_I$'),TeX('$\\theta_A$')),pos=3)
