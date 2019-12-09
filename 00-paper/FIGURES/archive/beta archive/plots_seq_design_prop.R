par(mfrow = c(1,1)) 
#c(bottom, left, top, right)
par(mar=c(5.1+5,4.1+1,2.1,2.1))
i<-1
plot(p.range,outer[i,,"eff.mon.initial"],type='l',ylim=c(0,1),lwd=2,col="blue",
     ylab="Probability",xlab="",
     main=label_main,
     axes=FALSE)
box()
#text(p.range,outer[i,,"eff.mon.initial"],labels=format(round(outer[i,,"eff.mon.initial"],digits=2),nsmall=2),pos=3)
#lines(p.range,outer[i,,"fut.mon.final"],lwd=2)
lines(p.range,outer[i,,"fut.mon.initial"],lwd=2,col="red")
#text(p.range,outer[i,,"fut.mon.initial"],labels=format(round(outer[i,,"fut.mon.initial"],digits=2),nsmall=2),pos=1)
#lines(p.range,outer[i,,"eff.mon.final"],lwd=2)

## calculate futility
lines(p.range,1-outer[i,,"fut.mon.initial"]-outer[i,,"eff.mon.initial"],lwd=2,col="darkgrey")
#lines(p.range,outer[i,,"inc.final"],lwd=2)
#text(p.range,outer[i,,"inc"],
#     labels=format(round(outer[i,,"inc"],digits=2),nsmall=2),pos=1)

axis(1,las=0,at=p.range,labels=format(p.range,nsmall=2))
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')
abline(v=c(p.skpt,p.enth),col='grey',lty='dashed')

# legend("right",text.width=0.05,seg.len=0.2,
#        inset = c(0, -0.2), bty = "n", x.intersp=0.5,
#        xjust=0, yjust=0,
#        legend=c("Stop Early for Efficacy","Stop Early for Futility","Inconclusive Findings w/ Full Dataset"),
#        col=c("blue","red","darkgrey"),
#        lwd=3, cex = 0.75, xpd = TRUE)
row<-2
for (j in 1:length(p.range)){
  mtext(text=paste0(format(round(outer[i,j,"eff.mon.initial"],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[j])
}
mtext(text="EFF",side=1,line=row+1,at=stretch)
row<-3
for (j in 1:length(p.range)){
  mtext(text=paste0(format(round(outer[i,j,"fut.mon.initial"],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[j])
}
mtext(text="FUT",side=1,line=row+1,at=stretch)
row<-4
for (j in 1:length(p.range)){
  mtext(text=paste0(format(round(1-outer[i,j,"fut.mon.initial"]-outer[i,j,"eff.mon.initial"],
                                 digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[j])
}
mtext(text="INC",side=1,line=row+1,at=stretch)
row<-5
for (j in 1:length(p.range)){
  mtext(text=paste0(format(round(outer[i,j,"ss.initial"],digits=1),nsmall=1),
                    " + ",
                    format(round(outer[i,j,"ss.final"]-
                                   outer[i,j,"ss.initial"],digits=1),nsmall=1),
                    " = ",
                    format(round(outer[i,j,"ss.final"],digits=1),nsmall=1)),
        side=1,line=row+1,at=p.range[j])
}
mtext(text="SS",side=1,line=row+1,at=stretch)
row<-6
for (j in 1:length(p.range)){
  mtext(text=paste0("(I) ",
                    format(round(outer[i,j,"post.mean.initial"],digits=3),nsmall=3),
                    " (F) ",
                    format(round(outer[i,j,"post.mean.final"],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[j])
}
mtext(text="CP",side=1,line=row+1,at=stretch)
row<-7
for (j in 1:length(p.range)){
  mtext(text=paste0("(I) ",
                    format(round(outer[i,j,"cov.initial"],digits=3),nsmall=3),
                    " (F) ",
                    format(round(outer[i,j,"cov.final"],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[j])
}
mtext(text="PM",side=1,line=row+1,at=stretch)

text(x = c(0.365,0.215,0.25), 
     y = c(0.85,0.85,0.115),
     labels = c("Stop Early for Efficacy","Stop Early for Futility",
                "Inconclusive w/ Full Dataset"))


# legend(x=0.325,y=0.6,
#        legend=c("Stop Early for Efficacy",
#                 "Stop Early for Futility",
#                 "Inconclusive with Full Data"),
#        lty=c('longdash','dotted','solid'),
#        cex=0.5,
#        box.lwd = 1,box.col = "black",bg = "white",pt.cex = 1,
#        text.width=20)