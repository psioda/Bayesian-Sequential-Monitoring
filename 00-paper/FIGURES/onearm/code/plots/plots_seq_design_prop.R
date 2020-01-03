# ### Sequential design properties ########################
# label_main=""
stretch<-p.skpt*0.95
# pdf(file=paste0("seq_dsn",prior_indicator_skpt[idx],prior_indicator_enth[idx],".pdf"),height=6,width=15)
# source("plots/plots_seq_design_prop.R")
# dev.off()

par(mfrow = c(1,1)) #c(bottom, left, top, right)
par(mar=c(5.1+5,4.1+1,2.1,2.1))
i<-1
plot(figure3$p.range,figure3$eff.mon.initial,type='l',ylim=c(0,1),lwd=2,col="blue",
     ylab="Probability",xlab="",
     main=label_main,
     axes=FALSE)
box()
#text(figure3$figure3$p.range,outer[i,,"eff.mon.initial"],labels=format(round(outer[i,,"eff.mon.initial"],digits=2),nsmall=2),pos=3)
#lines(figure3$figure3$p.range,outer[i,,"fut.mon.final"],lwd=2)
lines(figure3$p.range,figure3$fut.mon.initial,lwd=2,col="red")
#text(figure3$figure3$p.range,outer[i,,"fut.mon.initial"],labels=format(round(outer[i,,"fut.mon.initial"],digits=2),nsmall=2),pos=1)
#lines(figure3$p.range,outer[i,,"eff.mon.final"],lwd=2)

## calculate futility
lines(figure3$p.range,1-figure3$fut.mon.initial-figure3$eff.mon.initial,lwd=2,col="darkgrey")
#lines(figure3$p.range,outer[i,,"inc.final"],lwd=2)
#text(figure3$p.range,outer[i,,"inc"],
#     labels=format(round(outer[i,,"inc"],digits=2),nsmall=2),pos=1)

axis(1,las=0,at=figure3$p.range[seq(1,length(figure3$p.range),by=3)],
     labels=format(figure3$p.range[seq(1,length(figure3$p.range),by=3)],nsmall=2))
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
for (j in seq(1,length(figure3$p.range),by=3)){
  mtext(text=paste0(format(round(figure3$eff.mon.initial[j],digits=3),nsmall=3)),
        side=1,line=row+1,at=figure3$p.range[j])
}
mtext(text="EFF",side=1,line=row+1,at=stretch)
row<-3
for (j in seq(1,length(figure3$p.range),by=3)){
  mtext(text=paste0(format(round(figure3$fut.mon.initial[j],digits=3),nsmall=3)),
        side=1,line=row+1,at=figure3$p.range[j])
}
mtext(text="FUT",side=1,line=row+1,at=stretch)
row<-4
for (j in seq(1,length(figure3$p.range),by=3)){
  mtext(text=paste0(format(round(1-figure3$fut.mon.initial[j]-figure3$eff.mon.initial[j],
                                 digits=3),nsmall=3)),
        side=1,line=row+1,at=figure3$p.range[j])
}
mtext(text="INC",side=1,line=row+1,at=stretch)
row<-5
for (j in seq(1,length(figure3$p.range),by=3)){
  mtext(text=paste0(format(round(figure3$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(figure3$ss.final[j]-
                                 figure3$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(figure3$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row+1,at=figure3$p.range[j])
}
mtext(text="SS",side=1,line=row+1,at=stretch)
row<-6
for (j in seq(1,length(figure3$p.range),by=3)){
  mtext(text=paste0("(I) ",
                    format(round(figure3$post.mean.initial[j],digits=3),nsmall=3),
                    " (F) ",
                    format(round(figure3$post.mean.final[j],digits=3),nsmall=3)),
        side=1,line=row+1,at=figure3$p.range[j])
}
mtext(text="PM",side=1,line=row+1,at=stretch)
row<-7
for (j in seq(1,length(figure3$p.range),by=3)){
  mtext(text=paste0("(I) ",
                    format(round(figure3$cov.initial[j],digits=3),nsmall=3),
                    " (F) ",
                    format(round(figure3$cov.final[j],digits=3),nsmall=3)),
        side=1,line=row+1,at=figure3$p.range[j])
}
mtext(text="CP",side=1,line=row+1,at=stretch)

# text(x = c(0.365,0.215,0.25), 
#      y = c(0.85,0.85,0.115),
#      labels = c("Stop Early for Efficacy","Stop Early for Futility",
#                 "Inconclusive w/ Full Dataset"))


# legend(x=0.325,y=0.6,
#        legend=c("Stop Early for Efficacy",
#                 "Stop Early for Futility",
#                 "Inconclusive with Full Data"),
#        lty=c('longdash','dotted','solid'),
#        cex=0.5,
#        box.lwd = 1,box.col = "black",bg = "white",pt.cex = 1,
#        text.width=20)