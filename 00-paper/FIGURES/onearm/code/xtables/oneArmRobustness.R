rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

## load data from longleaf
Table1    <- read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1 <- merge(args_simulation,Table1,by.x="X",by.y="idx")

combined1 <- combined1[(combined1$model==1 & combined1$p.range %in% c(0.4, 0.535, 0.67)),]
combined1$comb <- NA
combined1$comb[(combined1$skpt_spike==0 & combined1$enth_flat==0)]=1
combined1$comb[(combined1$skpt_spike==1 & combined1$enth_flat==0)]=2
combined1$comb[(combined1$skpt_spike==0 & combined1$enth_flat==1)]=3
combined1$comb[(combined1$skpt_spike==1 & combined1$enth_flat==1)]=4

combined1 <- combined1[,c("p.range","comb","fut.mon.initial","fut.mon.final","ss.initial")]
library(tidyr)
spread(combined1, comb, c(fut.mon.initial,fut.mon.final,ss.initial))




#par(mar=c(5.1+5,4.1+1,2.1,2.1))
#par(mai=c(1.2,1.1,0,1)) #c(bottom, left, top, right)
label_main=""
stretch<-p.skpt*0.9
width.scale<-7

png('../../figureS2a.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)
## subset data before running plot
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==1 & combined1$enth_flat==0,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
dev.off()

png('../../figureS2b.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,res=300)
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==1 & combined1$enth_flat==1,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
dev.off()

png('../../figureS2c.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,res=300)
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==0 & combined1$enth_flat==1,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
dev.off()

png('../../figureS2d.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,res=300)
figure3<-combined1[combined1$model==1 & combined1$skpt_spike==0 & combined1$enth_flat==0,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=3),]
source("plots/plots_seq_design_prop.R")
dev.off()