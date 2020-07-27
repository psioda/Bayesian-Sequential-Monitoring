plot(log2(figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]]),
     figure4$eff.mon.initial[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
     type='l',
     ylim=c(0,0.05),
     lwd=2,
     lty='longdash',
     ylab="Probability",
     xlab="",
     main="",
     axes=FALSE)
# horizontal lines
abline(h=seq(0,0.05,by=0.01),col='grey')
box()

for (i in length(k1)){
points(log2(figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]]),
      figure4$eff.mon.initial[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
      lwd=2)
points(log2(figure4$freq.mntr[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]]),
      figure4$eff.mon.final[figure4$out.mean==k1[1] & figure4$enr.shape==k2[1]],
      lwd=2,pch=20)
}

text(log2(figure4$freq.mntr),figure4$eff.mon.initial,
     labels=format(round(figure4$eff.mon.initial,digits=3),nsmall=3),pos=1)


lines(log2(figure4$freq.mntr),figure4$eff.mon.initial)
lines(log2(figure4$freq.mntr),figure4$eff.mon.final)

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
axis(1,las=0,at=log2(figure4$freq.mntr[figure4$freq.mntr %in% mntr.pts]),
     labels=figure4$freq.mntr[figure4$freq.mntr %in% mntr.pts])

# # manual x-axis
# mtext(text=figure4$freq.mntr,side=1,
#       line=1,
#       at=figure4$freq.mntr[figure4$out.mean==k1 & figure4$enr.shape==k2])
# 

mtext(text="Monitor Freq",side=1,line=1,at=stretch,adj=0)
mtext(text="SS",side=1,line=2,at=stretch,adj=0)
mtext(text="Follow-up",side=1,line=3,at=stretch,adj=0)
mtext(text="SS Final",side=1,line=4,at=stretch,adj=0)

# sample size information
mtext(text=paste0(
   format(round(figure4$ss.initial[figure4$freq.mntr %in% mntr.pts],digits=1),nsmall=1)),
      side=1,
      line=2,
      at=log2(figure4$freq.mntr[figure4$freq.mntr %in% mntr.pts]))

mtext(text=paste0(format(round(figure4$ss.final[figure4$freq.mntr %in% mntr.pts]-
                                   figure4$ss.initial[figure4$freq.mntr %in% mntr.pts],digits=1),nsmall=1)),
    side=1,
    line=3,
    at=log2(figure4$freq.mntr[figure4$freq.mntr %in% mntr.pts]))

mtext(text=paste0(format(round(figure4$ss.final[figure4$freq.mntr %in% mntr.pts],digits=1),nsmall=1)),
    side=1,
    line=4,
    at=log2(figure4$freq.mntr[figure4$freq.mntr %in% mntr.pts]))

