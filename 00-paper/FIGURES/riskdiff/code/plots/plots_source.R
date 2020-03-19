#########################################
#### Figure 5, Risk Diff Prior Plots ####
#########################################
rm(list = ls())

library(pracma)

load(file = '../args_model.RData') # loads all model information include prior parameters

source("../code_functions.R")
source("../code_integrate.R")

width.scale <- 6
source("figure5.R")

png('../../../figure5a.png',
    width = 300*width.scale, 
    height = 300*width.scale,
    pointsize = 16,
    res = 300)
source("figure5a.R")
mtext("(B)",
      side = 2,
      line = 1,
      at = 1,
      las = 1)
text(x = 0.9,
     y = 0.7,
     as.expression(bquote(theta[1]==theta[0]+delta[S])))
text(x = 0.5,
     y = 0.8,
     as.expression(bquote(theta[1]==theta[0]+delta[E])))
dev.off()

png('../../../figure5b.png',
    width = 300*width.scale, 
    height = 300*width.scale,
    pointsize = 16,
    res = 300)
source("figure5b.R")
mtext("(D)",
      side=2,
      line=1,
      at=1,
      las=1)
text(x=0.9,
     y=0.7,
     as.expression(bquote(theta[1]==theta[0]+delta[S])))
text(x=0.5,
     y=0.8,
     as.expression(bquote(theta[1]==theta[0]+delta[E])))
dev.off()

png('../../../figure5c.png',
    width = 300*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)
source("figure5c.R")
mtext("(A)",side=2,line=1,at=2.1,las=1)
dev.off()

png('../../../figure5d.png',
    width = 300*width.scale, 
    height = 300*width.scale,
    pointsize = 16,
    res = 300)
source("figure5d.R")
mtext("(C)",
      side = 2,
      line = 1,
      at = 9.25,
      las = 1)
dev.off()

png('../../../figure5e.png',
    width = 300*width.scale, 
    height = 300*width.scale,
    pointsize = 16,
    res = 300)
source("figure5e.R")
mtext("(E)",
      side = 2,
      line = 1,
      at = 9.25,
      las = 1)
dev.off()


#############################################
#### Figure 6, Risk Diff Inference Plots ####
#############################################

rm(list = ls())
width.scale<-7
png('../../../figure6.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)
par(mar=c(6.1+1,4.1+1,2.1,2.1)) #c(bottom, left, top, right)

stretch <- 0.37

load(file = '../args_model.RData') # loads all model information include prior parameters
args_simulation <- read.csv(file = "../args_simulation.csv", header = TRUE, sep = ",")

Table1 <- read.csv(file = "../../output/Table1_merged.csv", header = TRUE, sep = ",")
combined1 <- merge(args_simulation, Table1, by.x = "X", by.y = "idx")
figure3 <- combined1

plot(figure3$p.IP[figure3$eff.mix.prob==0.25],
     figure3$eff.mon.initial[figure3$eff.mix.prob==0.25],
     type='l',
     ylim=c(0,1),
     lwd=1,
     ylab="Probability",
     xlab="",
     main="",
     xaxt="n")
box()

probs <- seq(0.25, 1, by = 0.25)
for (i in 1:length(probs)){
  row <- i + 2
  
  temp <- figure3[is.na(figure3$eff.mix.prob) == FALSE,]
  temp <- temp[temp$eff.mix.prob == probs[i],]
  
  lines(temp$p.IP,temp$eff.mon.initial)
  
  for (j in seq(1,length(temp$p.IP), by=3)){
    mtext(text=paste0(#format(round(temp$ss.initial[j],digits=1),nsmall=1),
                      #" + ",
                      #format(round(temp$ss.final[j]-
                      #               temp$ss.initial[j],digits=1),nsmall=1),
                      #" = ",
                      format(round(temp$ss.final[j],digits=1),nsmall=1)),
          side=1,line=row,at=temp$p.IP[j])
  }
  text(temp$p.IP[j], temp$eff.mon.initial[j], i + 1)
  mtext(text=paste0(i + 1),side=1,line=row,at=stretch,adj=0)
}

#### NOW FOR NA EFF.MON.PROB ####
for (i in 5){
  row <- i - 3
  
  temp <- figure3[is.na(figure3$eff.mix.prob) == TRUE,]
  
  lines(temp$p.IP,temp$eff.mon.initial)
  
  for (j in seq(1,length(temp$p.IP), by=3)){
    mtext(text=paste0(#format(round(temp$ss.initial[j],digits=1),nsmall=1),
      #" + ",
      #format(round(temp$ss.final[j]-
      #               temp$ss.initial[j],digits=1),nsmall=1),
      #" = ",
      format(round(temp$ss.final[j],digits=1),nsmall=1)),
      side=1,line=row,at=temp$p.IP[j])
  }
  mtext(text=paste0(1),side=1,line=row,at=stretch,adj=0)
  text(temp$p.IP[j], temp$eff.mon.initial[j], 1)
  
}

axis(1,las=0,at=temp$p.IP[seq(1,length(temp$p.IP),by=3)],
     labels=format(temp$p.IP[seq(1,length(temp$p.IP),by=3)],nsmall=2))

#axis(2,las=0,at=seq(0,1,by=0.1), labels=seq(0,1,by=0.1))

mtext(text=c(as.expression(bquote(theta))),side=1,line=1,at=stretch,adj=0)

legend('topleft',
       legend= c("1: Adaptive Weight Mixture",
                 "2: 75:25 Mixture (Mostly Enthusiastic)",
                 "3: 50:50 Mixture",
                 "4: 25:75 Mixture",
                 "5: 100:0 Mixture (Default, All Skeptical)"))

dev.off()
