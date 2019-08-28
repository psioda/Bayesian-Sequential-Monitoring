
## POWER CURVES
par(ask=FALSE)
par(mfrow = c(1,1))
k<-1
plot(p.range,outer.eff[k,],type='l',ylim=c(0,1),lwd=2,
     ylab="Probability",xlab="",
     main="Power Curves (Probability of Efficacy at Interim)\n Modifying Frequency of Monitoring",
     axes=FALSE)
for (k in 2:6){
  lines(p.range,outer.eff[k,],lwd=2)
}
box()
axis(1,las=0,at=p.range,labels=format(p.range,nsmall=2))
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')
abline(v=c(p.skpt,p.enth),col='grey',lty='dashed')
row<-1
for (column in 1:length(p.range)){
  mtext(text=paste0(format(round(outer.ss.initial[k,column],digits=1),nsmall=1),
                    " + ",
                    format(round(outer.ss.final[k,column]-
                                   outer.ss.initial[k,column],digits=1),nsmall=1),
                    " = ",
                    format(round(outer.ss.final[k,column],digits=1),nsmall=1)),
        side=1,line=row+1,at=p.range[column])
}
row<-2
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",
                    format(round(outer.post.mean.initial[k,column],digits=3),nsmall=3),
                    " (F) ",
                    format(round(outer.post.mean.final[k,column],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[column])
}
row<-3
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",
                    format(round(outer.cov.initial[k,column],digits=3),nsmall=3),
                    " (F) ",
                    format(round(outer.cov.final[k,column],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[column])
}
column<-1
#legend(x=0.325,y=0.6,legend=c("Skeptic at Interim","Skeptical at Final","Mixture Prior"),
#       lty=c('longdash','dotted','solid'),
#       cex=0.8,
#       box.lwd = 1,box.col = "black",bg = "white",pt.cex = 1)
mtext(text="SS",side=1,line=2,at=0.125)
mtext(text="PM",side=1,line=3,at=0.125)
mtext(text="CP",side=1,line=4,at=0.125)

## Plot of power curves
par(ask=TRUE)
plot(p.range,outer_trial_result_binary[1,],type='l',
     xlab="Response Probability",ylab="Probability of Rejection",main="Power Curves")
lines(p.range,outer_trial_result_binary[2,],type='l')
lines(p.range,outer_trial_result_binary[3,],type='l')
lines(p.range,outer_trial_result_binary[4,],type='l')
legend(x=p.enth*(2/3),y=0.8,
       title="# enrolled",
       legend=freq.mntr,
       col=freq.mntr,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=TRUE)
plot(p.range,outer.ss.initial[4,],type='l',
     xlab="Response Probability",ylab="Expected Sample Size",main="Sample Size",
     ylim=c(max.ss/2,max.ss*1.1))
axis(2,at=max.ss,labels="Max")
points(p.range,outer.ss.initial[3,],type='l',lty=2)
points(p.range,outer.ss.initial[2,],type='l',lty=3)
points(p.range,outer.ss.initial[1,],type='l',lty=4)
legend(x=p.skpt*(3/2),max.ss*(4/5),
       title="# enrolled",
       legend=freq.mntr,
       col=freq.mntr,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=TRUE)
plot(p.range,abs(outer.phat.initial[4,]-p.range),type='l',
     xlab="Response Probability",ylab="Bias",main="Bias",
     ylim=c(0,0.1))
points(p.range,abs(outer.phat.initial[3,]-p.range),type='l',lty=2)
points(p.range,abs(outer.phat.initial[2,]-p.range),type='l',lty=3)
points(p.range,abs(outer.phat.initial[1,]-p.range),type='l',lty=4)
legend(x=p.enth*(6/7),0.1,
       title="# enrolled",
       legend=freq.mntr,
       col=freq.mntr,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=FALSE)
