#######################################
#### Figure 1, One Arm Prior Plots ####
#######################################
dev.off()
rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)
source("code_functions.R")
source("args_model.R")
width.scale<-6
png(paste0(plot.dir,'figure1a.png'),width = 300*width.scale, height = 300*width.scale,pointsize=16,res=300)
prior.nc.skpt<-skpt_prior_default()
source("plots/plots_prior_skpt.R")
mtext("(A)",side=2,line=1,at=6,las=1)
dev.off()

rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)
source("code_functions.R")
source("args_model.R")
width.scale<-6
png(paste0(plot.dir,'figure1b.png'),width = 300*width.scale, height = 300*width.scale,pointsize=16,res=300)
prior.nc.skpt<-skpt_prior_custom(scale=1.15)
source("plots/plots_prior_skpt.R")
mtext("(B)",side=2,line=1,at=6,las=1)
dev.off()

rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)
source("code_functions.R")
source("args_model.R")
width.scale<-6
png(paste0(plot.dir,'figure1c.png'),width = 300*width.scale, height = 300*width.scale,pointsize=16,res=300)
prior.nc.enth<-enth_prior_default()
source("plots/plots_prior_enth.R")
mtext("(C)",side=2,line=1,at=3,las=1)
dev.off()

rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)
source("code_functions.R")
source("args_model.R")
width.scale<-6
png(paste0(plot.dir,'figure1d.png'),width = 300*width.scale, height = 300*width.scale,pointsize=16,res=300)
prior.nc.enth<-enth_prior_custom(scale=0.85)
source("plots/plots_prior_enth.R")
mtext("(D)",side=2,line=1,at=3,las=1)
dev.off()

########################################
#### Figure 2, One Arm Violin Plots ####
########################################
rm(list = ls())
plot.dir <-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/violin plots and others/"

root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)
source("code_functions.R")
source("args_model.R")
prior.nc.skpt<-skpt_prior_custom(scale=0.75)
prior.nc.enth<-enth_prior_default()
label.x<- (-15)
x.len<-1000
grid<-seq(0+1E-4,1-1E-4,length=x.len)
width.scale<-7
png(paste0(plot.dir,'figure2a.png'),width = 450*width.scale, height = 300*width.scale,pointsize=16,res=300)
par(mar=c(5.1+1,4.1,2.1,2.1)) #c(bottom, left, top, right)
# Efficacy example
miss<-c(0,2,3  ,4 ,0)
n<-   c(0,10,20,30,34)
y1<-  c(0,7,12 ,21,21)
y0=n-y1
spacing<-(seq(1:length(n))-1)*12
source("plots/plots_violin_v2.R")
legend("bottomright",
       c("Skeptical Posterior",
         "Enthuastic Posterior"),
       #bty="n",
       fill=c("darkgray","lightgray"),
       density=c(NA,NA),
       lty=c(NA,NA),
       pch=c(NA,NA),
       lwd=c(NA,NA),
       border=c('black','black')#,x.intersp=c(-.5,-.5,1,1)
)
mtext("(A)", side = 2, line = 3, at = 1, las = 1)
dev.off()

rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
plot.dir <-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/violin plots and others/"
setwd(root)
source("code_functions.R")
source("args_model.R")
prior.nc.skpt<-skpt_prior_custom(scale=0.75)
prior.nc.enth<-enth_prior_default()
label.x<- (-15)
x.len<-1000
grid<-seq(0+1E-4,1-1E-4,length=x.len)
width.scale<-7
png(paste0(plot.dir,'figure2b.png'),width = 450*width.scale, height = 300*width.scale,pointsize=16,res=300)
par(mar=c(5.1+1,4.1,2.1,2.1)) #c(bottom, left, top, right)
# Futility example
n<-   c(0,10,20,30,37)
y1<-  c(0,6 ,10 ,14 ,15)
miss<-c(0,2 ,4 ,7 ,0)
y0=n-y1
spacing<-(seq(1:length(n))-1)*12
source("plots/plots_violin_v2.R")
legend("topright",
       c("Skeptical Posterior",
         "Enthuastic Posterior"),
       #bty="n",
       fill=c("darkgray","lightgray"),
       density=c(NA,NA),
       lty=c(NA,NA),
       pch=c(NA,NA),
       lwd=c(NA,NA),
       border=c('black','black')#,x.intersp=c(-.5,-.5,1,1)
       )
mtext("(B)", side = 2, line = 3, at = 1, las = 1)
dev.off()


#########################################################
#### Figure 3a, One Arm Sequential Design Properties ####
#########################################################
dev.off()
par(mar=c(5.1,4.1,2.1,2.1)) #c(bottom, left, top, right)
rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

figure3<-combined1[combined1$model==1 & combined1$skpt_spike==1 & combined1$enth_flat==0,]
figure3_table<-figure3[seq(1,length(figure3$p.range),length=5),]

label_main=""
stretch<-p.skpt*0.9
width.scale<-7

png('../../figure3a.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)
source("plots/plots_seq_design_prop.R")
dev.off()

##############################################
#### Figure 3b, One Arm Evidence Decrease ####
##############################################
dev.off()
rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table2<-read.csv(file="../output/Table2_merged.csv",header=TRUE,sep=",")
combined2<-merge(args_simulation,Table2,by.x="X",by.y="idx")
Table3<-read.csv(file="../output/Table3_merged.csv",header=TRUE,sep=",")

# need to fix warning message
combined3<-merge(combined2,Table3,by.x="X",by.y="idx")

figure3b<-combined3[combined3$model==1 & combined3$skpt_spike==1 & combined3$enth_flat==0,]
figure3b_table<-figure3b[figure3b$p.range %in% c(figure3b$p.range[seq(1,length(figure3b$p.range),length=5)]),]
width.scale<-7


png('../../figure3b.png',width = 450*width.scale, height = 300*width.scale,pointsize=16,res=300)
source("plots/plots_decrease.R")
dev.off()

###############################
#### Figure 4, One Arm T1E ####
###############################
rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

## load data from longleaf
Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

## subset data before running plot
k1<-c(4)
k2<-c(1)
figure4<-combined1[combined1$model==2 & combined1$out.mean==k1[1] & combined1$enr.shape==k2[1],]
v<-floor(max.ss/figure4$freq.mntr)
#freq.mntr.lst<-length(v)-match(unique(v),rev(v))+1
#figure4<-figure4[figure4$freq.mntr %in% freq.mntr.lst,]

width.scale<-7

# run plot
#png('../../figure4.png',width = 450*width.scale, height = 300*width.scale,pointsize=16,res=300)
label_main=""
stretch<- (-1.15)

mntr.pts<-c(1,2,4,8,16,56,112)

source("plots/plots_t1e.R")
legend("right",c("Initial","Final"),lty=c(5,1))

#dev.off()


##################################################
### Figure S1, NEW One Arm Enrollment Schemes ####
##################################################
rm(list = ls())
root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")
k1<-c(4)
k2<-c(0.25)
width.scale<-7

## subset data before running plot, called figure4 for historical reasons
figure4<-combined1[combined1$model==2 & combined1$enr.shape==k2[1],]

# run plot
png('../../figureS1.png',width = 450*width.scale, height = 300*width.scale,pointsize=16,res=300)
label_main=""
stretch<- (-1.15)

mntr.pts<-c(1,2,4,8,16,56,112)
source("plots/plots_t1e.R")
legend("right",c("Initial","Final"),lty=c(5,1))

dev.off()

###############################################
#### Figure S1, One Arm Enrollment Schemes #### 
###############################################
# rm(list = ls())
# root<-"P:/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
# setwd(root)
# 
# source("args_model.R")
# args_simulation<-read.csv(file="args_simulation_S1.csv",header=TRUE,sep=",")[seq(141,188),]
# 
# Table1<-read.csv(file="../output/Table1_merged_S1.csv",header=TRUE,sep=",")
# combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")
# 
# figureS1<-combined1[combined1$model==3 & combined1$p.range==p.skpt,]
# 
# figureS1$id<-paste0(figureS1$enr.shape)
# 
# figureS1$eff<-format(round(figureS1$eff.mon.initial,digits=3),nsmall=3)
# figureS1$T1E<-format(round(figureS1$eff.mon.final,digits=3),nsmall=3)
# figureS1$ss<-format(round(figureS1$ss.final,digits=1),nsmall=1)
# figureS1$ongoing<-paste0(format(round((figureS1$ss.final-figureS1$ss.initial)/figureS1$ss.final*100,digits=1),nsmall=1),"\\%")
# 
# a1<-reshape(figureS1[,c("freq.mntr","eff","id")],
#             timevar = "freq.mntr",
#             idvar="id",
#             direction = "wide")
# names(a1) <- gsub("eff.", "", names(a1))
# a2<-reshape(figureS1[,c("freq.mntr","T1E","id")],
#             timevar = "freq.mntr",
#             idvar="id",
#             direction = "wide")
# names(a2) <- gsub("T1E.", "", names(a2))
# a3<-reshape(figureS1[,c("freq.mntr","ss","id")],
#             timevar = "freq.mntr",
#             idvar="id",
#             direction = "wide")
# names(a3) <- gsub("ss.", "", names(a3))
# a4<-reshape(figureS1[,c("freq.mntr","ongoing","id")],
#             timevar = "freq.mntr",
#             idvar="id",
#             direction = "wide")
# names(a4) <- gsub("ongoing.", "", names(a4))
# 
# all<-rbind(a1,a2,a3,a4)
# all <- all[ order(row.names(all)), ]
# write.table(t(all), "../../figureS1.csv", sep="&", col.names=T)
#######################################
#### Figure S2, One Arm Robustness #### 
#######################################
# rm(list = ls())
# root<-"P:/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
# setwd(root)
# 
# source("args_model.R")
# args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")
# 
# Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
# combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")
# 
# figureS2<-combined1[combined1$model==4,]
# 
# figureS2$id<-c(rep(1,5),rep(2,5),rep(3,5),rep(4,5))
# 
# figureS2$eff<-format(round(figureS2$eff.mon.initial,digits=3),nsmall=3)
# figureS2$fut<-format(round(figureS2$fut.mon.initial,digits=3),nsmall=3)
# figureS2$inc<-format(round(1-figureS2$eff.mon.initial-figureS2$fut.mon.initial,digits=3),nsmall=3)
# figureS2$ss<-paste0(format(round(figureS2$ss.initial,digits=1),nsmall=1),
#                     " + ",
#                     format(round(figureS2$ss.final-
#                                    figureS2$ss.initial,digits=1),nsmall=1),
#                     " = ",
#                     format(round(figureS2$ss.final,digits=1),nsmall=1))
# 
# a1<-reshape(figureS2[,c("p.range","eff","id")],
#             timevar = "p.range",
#             idvar="id",
#             direction = "wide")
# names(a1) <- gsub("eff.", "", names(a1))
# a2<-reshape(figureS2[,c("p.range","fut","id")],
#             timevar = "p.range",
#             idvar="id",
#             direction = "wide")
# names(a2) <- gsub("fut.", "", names(a2))
# a3<-reshape(figureS2[,c("p.range","inc","id")],
#             timevar = "p.range",
#             idvar="id",
#             direction = "wide")
# names(a3) <- gsub("inc.", "", names(a3))
# a4<-reshape(figureS2[,c("p.range","ss","id")],
#             timevar = "p.range",
#             idvar="id",
#             direction = "wide")
# names(a4) <- gsub("ss.", "", names(a4))
# 
# 
# all<-rbind(a1,a2,a3,a4)
# labels<-c("EFF","FUT","INC","SS")
# all$labels<-rep(c("EFF","FUT","INC","SS"),4)
# all<-all[,c("id","labels", "0.4","0.4675","0.535","0.6025","0.67")]
# all<-all[ order(row.names(all)), ]
# write.table(all, "../../figureS2.csv", sep="&", col.names=T, row.names=F)

############################################ 
#### Figure S2, NEW: ONE ARM ROBUSTNESS #### 
############################################ 
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
