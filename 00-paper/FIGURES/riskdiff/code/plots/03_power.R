#############################################
#### Figure 6, Risk Diff Inference Plots ####
#############################################

rm(list = ls())

width.scale <- 7
png('../../../figure6.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)
par(mar=c(6.1+1,4.1+1,2.1,2.1)) #c(bottom, left, top, right)
stretch <- 0.37 # to add x-axis table under graph

# set initial plotting area
plot(NULL,
     type = 'l',
     xlim = c(.39, .63),
     ylim = c(0,1),
     lwd  = 1,
     ylab = "Probability",
     xlab = "",
     main = "",
     xaxt = "n")
mtext(text=c(as.expression(bquote(theta))),side=1,line=1,at=stretch,adj=0)
legend('topleft',
       legend= c("1: Adaptive Weight Mixture (25% maximum)",
                 "2: 25% Skeptical (Mostly Enthusiastic)",
                 "3: 50% Skeptical",
                 "4: 75% Skeptical Mixture",
                 "5: 100% Skeptical (Default)"))


# plot all fixed weight priors
Table1        <- read.csv(file = "../../output/table1031920.csv", header = T)
#Table1_ref    <- read.csv(file = "../../output/output062120/Table1/1Table1.csv")
#names(Table1) <- names(Table1_ref)
args          <- read.csv(file = "../../output/args_simulation031920.csv", header = F)
args_ref      <- read.csv(file = "../args_simulation.csv")
names(args)   <- names(args_ref)
combined1     <- merge(args, Table1, by.x = "X", by.y = "idx")
figure3       <- combined1


## may need to change eff.mix.prob to eff.mix.prob.x
probs <- seq(0.25, 1, by = 0.25)
for (i in 1:length(probs)){
  row  <- i + 2
  temp <- figure3[is.na(figure3$eff.mix.prob) == FALSE,]
  temp <- temp[temp$eff.mix.prob == probs[i],]
  temp <- temp[temp$p.IP != 0.51,]
  lines(temp$p.IP,temp$eff.mon.initial)
  for (j in seq(1,length(temp$p.IP)#, by=3
  )){
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

# plot adaptive mixing prior
Table0  <- read.csv(file = "../Table0_merged062120.csv")
figure3 <- aggregate(x   = Table0,
                     by  = list(Table0$p.IP, Table0$p.PC, Table0$eff.mix.prob.x),
                     FUN = mean)
for (i in 5){
  row <- i - 3
  temp <- figure3[figure3$eff.mix.prob.x == 10, ]
  #temp <- figure3[is.na(figure3$eff.mix.prob.x) == TRUE,]
  lines(temp$p.IP,temp$eff.mon.initial)
  for (j in seq(1,length(temp$p.IP)#, by=3
  )){
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
axis(1,las=0,at=temp$p.IP[seq(1,length(temp$p.IP)#,by=3
)],
labels=format(temp$p.IP[seq(1,length(temp$p.IP)#,by=3
)],nsmall=2))
dev.off()
