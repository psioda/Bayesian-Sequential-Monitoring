#######################################
#### Figure 1, One Arm Prior Plots ####
#######################################
rm(list = ls())
root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("code_functions.R")
source("args_model.R")

png('../../figure1.png')
par(mfrow = c(2,2)) #c(bottom, left, top, right)
prior.nc.skpt<-skpt_prior_default()
source("plots/plots_prior_skpt.R")
prior.nc.skpt<-skpt_prior_custom(scale=1.15)
source("plots/plots_prior_skpt.R")
prior.nc.enth<-enth_prior_default()
source("plots/plots_prior_enth.R")
prior.nc.enth<-enth_prior_custom(scale=0.85)
source("plots/plots_prior_enth.R")
dev.off()

########################################
#### Figure 2, One Arm Violin Plots ####
########################################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)
library(latex2exp)
library(gridExtra)

root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("code_functions.R")
source("args_model.R")
prior.nc.skpt<-skpt_prior_custom(scale=1.15)
prior.nc.enth<-enth_prior_default()

# Efficacy example
n<-c(0,5,10,15,20)
y1<-c(0,4,8,12,16)
miss<-c(0,0,0,0,0)
y0=n-y1
stretch<-1.65

source("plots/plots_violin.R")
plot1 <- ggplot_gtable(ggplot_build(p))
plot1$layout$clip[plot1$layout$name == "panel"] <- "off"

png('../../figure2a.png')
par(mfrow = c(1,1)) #c(bottom, left, top, right)
grid.arrange(plot1)
dev.off()

# Futility example
n<-c(0,10,20,30,40)
y1<-c(0,2,4,6,8)
miss<-c(0,0,0,0,0)
y0=n-y1
stretch<-1.65 # "Sample Size" labels xlim
source("plots/plots_violin.R")
plot2 <- ggplot_gtable(ggplot_build(p))
plot2$layout$clip[plot2$layout$name == "panel"] <- "off"

png('../../figure2b.png')
par(mfrow = c(1,1)) #c(bottom, left, top, right)
grid.arrange(plot2)
dev.off()

# plot side-by-side
#grid.arrange(plot1,plot2,ncol=2)

#########################################################
#### Figure 3a, One Arm Sequential Design Properties ####
#########################################################
rm(list = ls())
root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

figure3<-combined1[combined1$model==1,]
label_main=""
png('../../figure3a.png',height=(480*2)*.45,width=(480*5)*.45)
source("plots/plots_seq_design_prop.R")
dev.off()

##############################################
#### Figure 3b, One Arm Evidence Decrease ####
##############################################
rm(list = ls())
root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table2<-read.csv(file="../output/Table2_merged.csv",header=TRUE,sep=",")
combined2<-merge(args_simulation,Table2,by.x="X",by.y="idx")
Table3<-read.csv(file="../output/Table3_merged.csv",header=TRUE,sep=",")
combined3<-merge(combined2,Table3,by.x="X",by.y="idx")

png('../../figure3b.png',height=(480*2)*.25,width=(480*5)*.25)
source("plots/plots_decrease.R")
dev.off()

###############################
#### Figure 4, One Arm T1E ####
###############################
rm(list = ls())
root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

figure4<-combined1[combined1$model==2,]
label_main=""
png('../../figure4.png',height=480,width=480)
source("plots/plots_t1e.R")
dev.off()

###############################################
#### Figure S1, One Arm Enrollment Schemes #### 
###############################################
rm(list = ls())
root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation_S1.csv",header=TRUE,sep=",")[seq(141,188),]

Table1<-read.csv(file="../output/Table1_merged_S1.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

figureS1<-combined1[combined1$model==3 & combined1$p.range==p.skpt,]

figureS1$id<-paste0(figureS1$enr.shape)

figureS1$eff<-format(round(figureS1$eff.mon.initial,digits=3),nsmall=3)
figureS1$T1E<-format(round(figureS1$eff.mon.final,digits=3),nsmall=3)
figureS1$ss<-format(round(figureS1$ss.final,digits=1),nsmall=1)
figureS1$ongoing<-paste0(format(round((figureS1$ss.final-figureS1$ss.initial)/figureS1$ss.final*100,digits=1),nsmall=1),"\\%")

a1<-reshape(figureS1[,c("freq.mntr","eff","id")],
            timevar = "freq.mntr",
            idvar="id",
            direction = "wide")
names(a1) <- gsub("eff.", "", names(a1))
a2<-reshape(figureS1[,c("freq.mntr","T1E","id")],
            timevar = "freq.mntr",
            idvar="id",
            direction = "wide")
names(a2) <- gsub("T1E.", "", names(a2))
a3<-reshape(figureS1[,c("freq.mntr","ss","id")],
            timevar = "freq.mntr",
            idvar="id",
            direction = "wide")
names(a3) <- gsub("ss.", "", names(a3))
a4<-reshape(figureS1[,c("freq.mntr","ongoing","id")],
            timevar = "freq.mntr",
            idvar="id",
            direction = "wide")
names(a4) <- gsub("ongoing.", "", names(a4))

all<-rbind(a1,a2,a3,a4)
all <- all[ order(row.names(all)), ]
write.table(t(all), "../../figureS1.csv", sep="&", col.names=T)
#######################################
#### Figure S2, One Arm Robustness #### 
#######################################
rm(list = ls())
root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

Table1<-read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1<-merge(args_simulation,Table1,by.x="X",by.y="idx")

figureS2<-combined1[combined1$model==4,]

figureS2$id<-c(rep(1,5),rep(2,5),rep(3,5),rep(4,5))

figureS2$eff<-format(round(figureS2$eff.mon.initial,digits=3),nsmall=3)
figureS2$fut<-format(round(figureS2$fut.mon.initial,digits=3),nsmall=3)
figureS2$inc<-format(round(1-figureS2$eff.mon.initial-figureS2$fut.mon.initial,digits=3),nsmall=3)
figureS2$ss<-paste0(format(round(figureS2$ss.initial,digits=1),nsmall=1),
                    " + ",
                    format(round(figureS2$ss.final-
                                   figureS2$ss.initial,digits=1),nsmall=1),
                    " = ",
                    format(round(figureS2$ss.final,digits=1),nsmall=1))

a1<-reshape(figureS2[,c("p.range","eff","id")],
            timevar = "p.range",
            idvar="id",
            direction = "wide")
names(a1) <- gsub("eff.", "", names(a1))
a2<-reshape(figureS2[,c("p.range","fut","id")],
            timevar = "p.range",
            idvar="id",
            direction = "wide")
names(a2) <- gsub("fut.", "", names(a2))
a3<-reshape(figureS2[,c("p.range","inc","id")],
            timevar = "p.range",
            idvar="id",
            direction = "wide")
names(a3) <- gsub("inc.", "", names(a3))
a4<-reshape(figureS2[,c("p.range","ss","id")],
            timevar = "p.range",
            idvar="id",
            direction = "wide")
names(a4) <- gsub("ss.", "", names(a4))


all<-rbind(a1,a2,a3,a4)
labels<-c("EFF","FUT","INC","SS")
all$labels<-rep(c("EFF","FUT","INC","SS"),4)
all<-all[,c("id","labels", "0.4","0.4675","0.535","0.6025","0.67")]
all<-all[ order(row.names(all)), ]
write.table(all, "../../figureS2.csv", sep="&", col.names=T, row.names=F)
