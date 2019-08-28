#### Type I error graphs ####
par(mar=c(5.1+4,4.1,4.1,2.1+1))
par(mfrow = c(1, 1)) 
k1<-4
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff"][out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='4 month follow-up, 2 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
mtext(text="SS",side=1,line=3,at=-.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff"][out.mean==k1 & enr.shape==k2],
     labels=format(round(outer[,,"eff"][out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer[,,"eff.final"][out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff.final"][out.mean==k1 & enr.shape==k2],
     labels=format(round(outer[,,"eff.final"][out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer[,,"ss.initial"][out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer[,,"ss.final"][out.mean==k1 & enr.shape==k2][column]-
                   outer[,,"ss.initial"][out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer[,,"ss.final"][out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=9-column,
    at=1/freq.mntr[column])
}
k1<-8
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='8 month follow-up, 2 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff.final[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff.final[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                   outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=9-column,
    at=1/freq.mntr[column])
}
par(mar=c(5.1+4,4.1,4.1,2.1+1))
k1<-4
k2<-0.25
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='4 month follow-up, 8 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff.final[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff.final[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
row<-1
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                   outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=9-column,
    at=1/freq.mntr[column])
}
k1<-8
k2<-0.25
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='8 month follow-up, 8 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff.final[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff.final[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                   outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=9-column,
    at=1/freq.mntr[column])
}
