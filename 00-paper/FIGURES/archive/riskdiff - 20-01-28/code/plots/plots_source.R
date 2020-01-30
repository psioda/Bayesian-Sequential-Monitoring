#########################################
#### Figure 5, Risk Diff Prior Plots ####
#########################################
rm(list = ls())
width.scale<-9
source("figure5.R")
source("figure5a.R")
mtext("(A)",side=2,line=3,at=1,las=1)
text(x=0.9,y=0.8,
     as.expression(bquote(theta[1]==theta[0]+delta[S])))
text(x=0.6,y=0.8,
     as.expression(bquote(theta[1]==theta[0]+delta[E])))
dev.off()


rm(list = ls())
width.scale<-9
source("figure5.R")
source("figure5b.R")
mtext("(B)",side=2,line=3,at=1,las=1)
text(x=0.9,y=0.8,
     as.expression(bquote(theta[1]==theta[0]+delta[S])))
text(x=0.6,y=0.8,
     as.expression(bquote(theta[1]==theta[0]+delta[E])))
dev.off()

#############################################
#### Figure 6, Risk Diff Inference Plots ####
#############################################
rm(list = ls())
root<-"P:/Git/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")
figure3<-combined1

par(mar=c(5.1,4.1,2.1,2.1))
par(mfrow = c(1,1))

png('../../figure6.png')
par(mfrow = c(1,1))

plot(figure3$p.IP,figure3$inference0,type='l',
     ylim=c(0,1),
     ylab="",
     xlab="",
     main="",
     axes=FALSE)
abline(h=seq(0,1,by=0.1),col='grey')
box()
lines(figure3$p.IP,figure3$inference0)
lines(figure3$p.IP,figure3$inference25)
lines(figure3$p.IP,figure3$inference50)
lines(figure3$p.IP,figure3$inference75)
lines(figure3$p.IP,figure3$inference100)
axis(1,las=0,at=figure3$p.IP[seq(1,length(figure3$p.IP),by=3)],
     labels=format(figure3$p.IP[seq(1,length(figure3$p.IP),by=3)],nsmall=2))
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))

dev.off()
