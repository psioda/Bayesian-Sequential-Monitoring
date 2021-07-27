rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

## load data from longleaf
Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

#par(mar=c(5.1+5,4.1+1,2.1,2.1))
#par(mai=c(1.2,1.1,0,1)) #c(bottom, left, top, right)
label_main=""
stretch<-p.skpt*0.9
width.scale<-7

png('plots/power-merged.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)

plot(NULL,
     xlim = c(mu0.skpt, mu0.enth),
     ylim = c(0, 1),
     ylab="Probability",
     xlab="",
     main=label_main,
     axes=FALSE)
box()
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')
axis(1,las=0,at=seq(mu0.skpt, mu0.enth, length = 3),
     labels=format(seq(mu0.skpt, mu0.enth, length = 3),nsmall=2))
mtext(bquote(theta),side=1,line=1,at=stretch)


## subset data before running plot
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==1 & combined1$enth_flat==0,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
figure3$eff.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$fut.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.final[figure3$p.range %in% c(0.4, 0.535, 0.67)]

source("plots/plots_seq_design_prop.R")

figure3<-combined1[combined1$model==1 & combined1$skpt_spike==1 & combined1$enth_flat==1,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
figure3$eff.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$fut.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.final[figure3$p.range %in% c(0.4, 0.535, 0.67)]
source("plots/plots_seq_design_prop.R")

figure3<-combined1[combined1$model==1 & combined1$skpt_spike==0 & combined1$enth_flat==1,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
figure3$eff.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$fut.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.final[figure3$p.range %in% c(0.4, 0.535, 0.67)]
source("plots/plots_seq_design_prop.R")

figure3<-combined1[combined1$model==1 & combined1$skpt_spike==0 & combined1$enth_flat==0,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
figure3$eff.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$fut.mon.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.initial[figure3$p.range %in% c(0.4, 0.535, 0.67)]
figure3$ss.final[figure3$p.range %in% c(0.4, 0.535, 0.67)]
source("plots/plots_seq_design_prop.R")

legend("right",#text.width=0.05,
       legend=c("Stop Early for Efficacy",
                "Stop Early for Futility"),
       lty=c('solid','longdash'))

dev.off()

png('plots/figureS2b.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,res=300)
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==1 & combined1$enth_flat==1,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
mtext("(B)", side = 2, line = 3, at = 1, las = 1)
dev.off()

png('plots/figureS2c.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,res=300)
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==0 & combined1$enth_flat==1,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
mtext("(C)", side = 2, line = 3, at = 1, las = 1)
dev.off()

png('plots/figureS2d.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,res=300)
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==0 & combined1$enth_flat==0,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
mtext("(D)", side = 2, line = 3, at = 1, las = 1)
dev.off()
