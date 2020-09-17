plot(figure3b$p.range,figure3b$X0.5,
     type='l',
     ylim=c(0.88,0.975),
     ylab="",
     xlab="",
     main="",
     axes=FALSE)
box()

axis(2,las=2,at=seq(0.88,1,by=0.01),labels=format(seq(0.88,1,by=0.01),nsmall=1))
title(ylab="Final Efficacy Criteria", line=3, #cex.lab=1.2
)
abline(h=seq(0.88,1,by=0.01),col='grey')

lines(figure3b$p.range,figure3b$X0.1)
lines(figure3b$p.range,figure3b$X0.25)
lines(figure3b$p.range,figure3b$X0.75)
lines(figure3b$p.range,figure3b$X0.9)
lines(c(p.skpt,p.enth),rep(0.975,2))

polygon(x=c(figure3b$p.range,rev(figure3b$p.range),figure3b$p.range[1]),
        y=c(figure3b$X0.1,rev(figure3b$X0.25),figure3b$X0.1[1]),
        col='lightgrey')

polygon(x=c(figure3b$p.range,rev(figure3b$p.range),figure3b$p.range[1]),
        y=c(figure3b$X0.9,rev(figure3b$X0.75),figure3b$X0.9[1]),
        col='lightgrey')

polygon(x=c(figure3b$p.range,rev(figure3b$p.range),figure3b$p.range[1]),
        y=c(figure3b$X0.75,rev(figure3b$X0.25),figure3b$X0.75[1]),
        col='darkgrey')

polygon(x=c(figure3b$p.range,rev(figure3b$p.range),figure3b$p.range[1]),
        y=c(figure3b$X0.9,rep(0.975,length(figure3b$p.range)),figure3b$X0.9[1]),
        col='black',density=15,angle=45,lty=1)

lines(figure3b$p.range,figure3b$X0.5)


## SUBSET DATA ##
axis(1,las=0,at=figure3b_table$p.range[seq(1,length(figure3b_table$p.range))],
     labels=format(figure3b_table$p.range[seq(1,length(figure3b_table$p.range))],nsmall=2))

text(figure3b_table$p.range,
     figure3b_table$X0.5,
     labels=format(round(figure3b_table$X0.5,digits=3),nsmall=3),
     pos=1)

points(figure3b_table$p.range,figure3b_table$X0.5,pch=20)

for (j in seq(1,length(figure3b_table$p.range))){
  #points(figure3b_table$p.range[j],figure3b_table$X0.1[j],pch=20)
  #text(figure3b_table$p.range[j]-0.005,figure3b_table$X0.1[j],
  #     labels=paste0(format(round(figure3b_table$X0.1[j],digits=3),nsmall=3)))
  
  #points(figure3b_table$p.range[j],figure3b_table$X0.25[j],pch=20)
  #text(figure3b_table$p.range[j]-0.005,figure3b_table$X0.25[j],
  #     labels=paste0(format(round(figure3b_table$X0.25[j],digits=3),nsmall=3)))
  
  #points(figure3b_table$p.range[j],figure3b_table$X0.5[j],pch=20)
  #text(figure3b_table$p.range[j]-0.005,figure3b_table$X0.5[j],
  #     labels=paste0(format(round(figure3b_table$X0.5[j],digits=3),nsmall=3)))
  
  #lines(x=rep(figure3b_table$p.range[j],2),c(figure3b_table$X0.1[j],0.975))
  
  mtext(text=paste0(format(round(figure3b_table$conditional[j],digits=3),nsmall=3)),
        side=1,
        line=2,
        at=figure3b_table$p.range[j])
}

lines(x=c(p.skpt,p.enth),rep(0.975,2))

#for (j in seq(1,length(figure3$p.range),by=3)){
#  mtext(text=paste0(format(round(figure3$eff.mon.initial[j],digits=3),nsmall=3)),
#        side=1,line=row+1,at=figure3$p.range[j])
#}
stretch<-0.37
mtext(bquote(theta),side=1,line=1,at=stretch)
mtext(text="% AGR",side=1,line=2,at=stretch)
leg <-legend("bottomright",
       c("Median","P10-P90","P25-75","P90-P100"), 
       fill=c("black","lightgray","darkgray","black"),
       density=c(0,NA,NA,10),
       pch=c(20,NA,NA,NA),
       lty=c(1,NA,NA,NA),
       border=c(NA,'black','black','black'),
       x.intersp=c(0,0,0,0))

# 
# legend("topright",
#        c("Skeptical Prior Posterior",
#          "Enthuastic Prior Posterior",
#          "Mixture Prior 95% Credible Interval",
#          "Mixture Prior Posterior Mean"),
#        #bty="n",
#        fill=c("darkgray","lightgray","black","black"),
#        density=c(NA,NA,0,0),
#        lty=c(NA,NA,1,1),
#        pch=c(NA,NA,3,20),
#        lwd=c(NA,NA,1,1),
#        border=c('black','black',NA,NA),
#        x.intersp=c(-.5,-.5,1,1))