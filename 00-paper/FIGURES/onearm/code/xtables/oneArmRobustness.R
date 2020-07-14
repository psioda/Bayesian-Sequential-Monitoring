rm(list = ls())

library(stringi)
library(xtable)

root<-"/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
setwd(root)

source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

## load data from longleaf
Table1    <- read.csv(file="../output/Table1_merged.csv",header=TRUE,sep=",")
combined1 <- merge(args_simulation,Table1,by.x="X",by.y="idx")

combined1      <- combined1[(combined1$model==1 & combined1$p.range %in% c(0.4, 0.535, 0.67)),]
combined1$comb <- NA
combined1$comb[(combined1$skpt_spike==0 & combined1$enth_flat==0)] = 1
combined1$comb[(combined1$skpt_spike==1 & combined1$enth_flat==0)] = 2
combined1$comb[(combined1$skpt_spike==0 & combined1$enth_flat==1)] = 3
combined1$comb[(combined1$skpt_spike==1 & combined1$enth_flat==1)] = 4

combined1 <- combined1[,c("p.range","comb","fut.mon.initial","fut.mon.final","ss.initial")]

M1        <- data.frame(t(combined1[combined1$p.range==0.4,]))
names(M1) <- M1["comb",]
M1$prob   <- M1[1,1]
M1        <- M1[!row.names(M1) %in% c("p.range","comb"),]

M2        <- data.frame(t(combined1[combined1$p.range==0.535,]))
names(M2) <- M2["comb",]
M2$prob   <- M2[1,1]
M2        <- M2[!row.names(M2) %in% c("p.range","comb"),]

M3        <- data.frame(t(combined1[combined1$p.range==0.67,]))
names(M3) <- M3["comb",]
M3$prob   <- M3[1,1]
M3        <- M3[!row.names(M3) %in% c("p.range","comb"),]

merged     <- rbind(M1,M2,M3)
merged$col <- row.names(merged)
merged     <- merged[,c("prob","col","1","2","3","4")]

for (i in 1:nrow(merged)){
    if (i%%nrow(M1)!=1) {merged$prob[i]=""}
    if (gsub('[[:digit:]]+', '', merged$col[i])=="fut.mon.initial") {merged$col[i]="FUT INITIAL"}
    if (gsub('[[:digit:]]+', '', merged$col[i])=="fut.mon.final") {merged$col[i]="FUT FINAL"}
    if (gsub('[[:digit:]]+', '', merged$col[i])=="ss.initial") {merged$col[i]="SS"}
}
print(xtable(merged, digits = c(rep(0,3),rep(4,4)), type = "latex"), 
      include.rownames=FALSE,
      include.colnames = FALSE,
      caption.placement = "top",
      file = paste0("../../4comparison.tex"),
      sanitize.text.function=function(x){x})

print(xtable(table, digits = c(rep(0,2),rep(4,4)), type = "latex", 
             caption = paste0("Operating characteristics for Design ", title_string,  " by historical sample size $n_0$ and power prior weight $a_0$.")), 
      include.rownames=FALSE,
      include.colnames = FALSE,
      file = paste0("output/xtables/sim", sim.vec[i], "-const", gsub("[.]","",const.vec[j]), ".tex"),       
      caption.placement = "top",
      sanitize.text.function=function(x){x}, 
      add.to.row=addtorow)
