#########################################
#### Figure 5, Risk Diff Prior Plots ####
#########################################
rm(list = ls())

library(pracma)
source("../args_model.R")
source("../code_functions.R")
source("../code_integrate.R")
source("../code_fcn_prior_placebo.R")
source("../code_skpt_tail_area.R")
source("../code_enth_tail_area.R")

fcn_prior_placebo()
skpt_tail_area()
enth_tail_area()
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


################## FIGURE 6 ###########################

rm(list = ls())
width.scale<-7
png('../../figure6.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)
par(mar=c(5.1+1,4.1+1,2.1,2.1)) #c(bottom, left, top, right)

stretch <- 0.365
#############################################
#### Figure 6, Risk Diff Inference Plots ####
#############################################
load(file = '../args_model.RData') # loads all model information include prior parameters
args_simulation <- read.csv(file = "../args_simulation.csv", header = TRUE, sep = ",")

Table1 <- read.csv(file = "../../output/Table1_merged.csv", header = TRUE, sep = ",")
combined1 <- merge(args_simulation, Table1, by.x = "X", by.y = "idx")
figure3 <- combined1

plot(figure3$p.IP[figure3$eff.mix.prob==0.25],
     figure3$eff.mon.initial[figure3$eff.mix.prob==0.25],
     type='l',
     ylim=c(0,0.5),
     lwd=1,
     ylab="Probability",
     xlab="",
     main="",
     xaxt="n")
box()



probs <- seq(0.25, 1, by = 0.25)
for (i in 1:length(probs)){
  row <- i + 1
  temp <- figure3[figure3$eff.mix.prob==probs[i],]
  
  lines(temp$p.IP,temp$eff.mon.initial)
  
  for (j in seq(1,length(temp$p.IP), by=3)){
    mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                      " + ",
                      format(round(temp$ss.final[j]-
                                     temp$ss.initial[j],digits=1),nsmall=1),
                      " = ",
                      format(round(temp$ss.final[j],digits=1),nsmall=1)),
          side=1,line=row,at=temp$p.IP[j])
  }
  mtext(text=c(as.expression(bquote(omega == .(probs[i])))),side=1,line=row,at=stretch,adj=0)
  
}

axis(1,las=0,at=temp$p.IP[seq(1,length(temp$p.IP),by=3)],
     labels=format(temp$p.IP[seq(1,length(temp$p.IP),by=3)],nsmall=2))

legend('topleft',
       legend= c(as.expression(bquote(omega == 0.25)),
                 as.expression(bquote(omega == 0.5)),
                 as.expression(bquote(omega == 0.75)),
                 as.expression(bquote(omega == 1))))


dev.off()

