#### Type I error graphs ####
#c(bottom, left, top, right)
par(mar=c(5.1+4,4.1+2,2.1,2.1))
par(mfrow = c(1, 1)) 
k1<-c(4)
k2<-c(1)
plot(rev(figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]]),
     figure4$eff.mon.initial[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
     type='l',
     ylim=c(0,0.05),
     lwd=2,
     lty='longdash',
     ylab="Type 1 Error Rate",
     xlab="",
     main="",
     axes=FALSE)
# horizontal lines
abline(h=seq(0,0.05,by=0.01),col='grey')
box()

for (i in length(k1)){
lines(rev(figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]]),
      figure4$eff.mon.final[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
      lwd=2)
}

# text(1/figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
#      figure4$eff.mon.initial[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
#      labels=format(round(figure4$eff.mon.initial[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],digits=3),
#                    nsmall=3),pos=3)
# text(1/figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
#      figure4$eff.mon.final[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
#      labels=format(round(figure4$eff.mon.final[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],digits=3),
#                    nsmall=3),pos=3)

# yaxixs
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
# x axis
axis(1,las=0,
     at=c(1,seq(10,100,by=10),max.ss),
     labels=rev(c(1,seq(10,100,by=10),max.ss)))

# # manual x-axis
# mtext(text=figure4$freq.mntr,side=1,
#       line=1,
#       at=figure4$freq.mntr[figure4$out.mean==k1 & figure4$enr.shape==k2])
# 

mtext(text="Monitor Freq",side=1,line=1,at=-15)
mtext(text="SS Interim",side=1,line=2,at=-15)
mtext(text="Follow-up",side=1,line=3,at=-15)
mtext(text="SS Final",side=1,line=4,at=-15)

# sample size information
mtext(text=paste0(format(round(figure4$ss.initial[c(1,seq(10,100,by=10),max.ss)],digits=1),nsmall=1)),
      side=1,
      line=2,
      at=rev(figure4$freq.mntr[c(1,seq(10,100,by=10),max.ss)]))

mtext(text=paste0(format(round(figure4$ss.final[c(1,seq(10,100,by=10),max.ss)]-
                               figure4$ss.initial[c(1,seq(10,100,by=10),max.ss)],digits=1),nsmall=1)),
      side=1,
      line=3,
      at=rev(figure4$freq.mntr[c(1,seq(10,100,by=10),max.ss)]))

mtext(text=paste0(format(round(figure4$ss.final[c(1,seq(10,100,by=10),max.ss)],digits=1),nsmall=1)),
      side=1,
      line=4,
      at=rev(figure4$freq.mntr[c(1,seq(10,100,by=10),max.ss)]))


