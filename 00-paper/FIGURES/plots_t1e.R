#### Type I error graphs ####
#c(bottom, left, top, right)
par(mar=c(5.1+4,4.1+2,2.1,2.1))
par(mfrow = c(1, 1)) 
k1<-4
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff.mon.initial"][out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main="",
     axes=FALSE)
# horizontal lines
abline(h=seq(0,0.15,by=0.01),col='grey')

text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff.mon.initial"][out.mean==k1 & enr.shape==k2],
     labels=format(round(outer[,,"eff.mon.initial"][out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer[,,"eff.mon.final"][out.mean==k1 & enr.shape==k2],lwd=2)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
outer[,,"eff.mon.initial"][out.mean==k1 & enr.shape==k2],lty='longdash',lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff.mon.final"][out.mean==k1 & enr.shape==k2],
     labels=format(round(outer[,,"eff.mon.final"][out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
box()


# yaxixs
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
# x axis
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=c("","","","","",""))

# manual x-axis
mtext(text=freq.mntr,side=1,line=1,at=1/freq.mntr[out.mean==k1 & enr.shape==k2])


# sample size information
mtext(text=paste0(format(round(
outer[,,"ss.initial"][out.mean==k1 & enr.shape==k2],
digits=1),nsmall=1)),
side=1,line=2,at=1/freq.mntr)
mtext(text=paste0(format(round(
outer[,,"ss.final"][out.mean==k1 & enr.shape==k2]-
outer[,,"ss.initial"][out.mean==k1 & enr.shape==k2],
digits=1),nsmall=1)),
side=1,line=3,at=1/freq.mntr)
mtext(text=paste0(format(round(
outer[,,"ss.final"][out.mean==k1 & enr.shape==k2],
digits=1),nsmall=1)),
side=1,line=4,at=1/freq.mntr)



mtext(text="Monitor Freq",side=1,line=1,at=-0.05)
mtext(text="SS Interim",side=1,line=2,at=-0.05)
mtext(text="Follow-up",side=1,line=3,at=-0.05)
mtext(text="SS Final",side=1,line=4,at=-0.05)

