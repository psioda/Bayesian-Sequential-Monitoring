figure3b<-combined3[combined3$model==1,]
figure3b<-figure3b[figure3b$p.range %in% c(figure3b$p.range[seq(1,length(figure3b$p.range),by=3)]),]

par(mfrow = c(1,1)) #c(bottom, left, top, right)
par(mar=c(5.1+0,4.1+0,2.1,2.1))

plot(figure3b$p.range,figure3b$X0.5,
     pch=20,
     ylim=c(0.9,0.975),
     ylab="",
     xlab="",
     main="",
     axes=FALSE)
box()

# xaxis
axis(1,las=0,at=figure3b$p.range[seq(1,length(figure3b$p.range))],
     labels=format(figure3b$p.range[seq(1,length(figure3b$p.range))],nsmall=2))
# yaxis
axis(2,las=2,at=seq(0.9,1,by=0.01),
     labels=format(seq(0.9,1,by=0.01),nsmall=1))
abline(h=seq(0.9,1,by=0.01),col='grey')

for (j in seq(1,length(figure3b$p.range))){
  points(figure3b$p.range[j],figure3b$X0.1[j],pch=20)
  #text(figure3b$p.range[j]-0.005,figure3b$X0.1[j],
  #     labels=paste0(format(round(figure3b$X0.1[j],digits=3),nsmall=3)))
  
  points(figure3b$p.range[j],figure3b$X0.25[j],pch=20)
  #text(figure3b$p.range[j]-0.005,figure3b$X0.25[j],
  #     labels=paste0(format(round(figure3b$X0.25[j],digits=3),nsmall=3)))
  
  points(figure3b$p.range[j],figure3b$X0.5[j],pch=20)
  #text(figure3b$p.range[j]-0.005,figure3b$X0.5[j],
  #     labels=paste0(format(round(figure3b$X0.5[j],digits=3),nsmall=3)))
  
  lines(x=rep(figure3b$p.range[j],2),c(figure3b$X0.1[j],0.975))
  
  mtext(text=paste0(format(round(figure3b$conditional[j],digits=3),nsmall=3)),
        side=1,
        line=3,
        at=figure3b$p.range[j])
}

lines(x=c(p.skpt,p.enth),rep(0.975,2))

#for (j in seq(1,length(figure3$p.range),by=3)){
#  mtext(text=paste0(format(round(figure3$eff.mon.initial[j],digits=3),nsmall=3)),
#        side=1,line=row+1,at=figure3$p.range[j])
#}
stretch<-0.36
mtext(text="% AGR",side=1,line=3,at=stretch)
