rm(list = ls())

#for (idx in 1:170){
if (.Platform$OS.type == "unix"){
  args<-commandArgs(trailingOnly = TRUE)
  idx<-as.numeric(args[1])}

if (.Platform$OS.type == "windows") {
  root<-"P:/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/code"
  setwd(root)
  idx<-170}

source("code_functions.R")
source("args_model.R")
args_simulation<-read.csv(file="args_simulation.csv",header=TRUE,sep=",")

# using idx to identify model
p.range<-args_simulation$p.range[idx]
freq.mntr<-args_simulation$freq.mntr[idx]
enr.shape<-args_simulation$enr.shape[idx]
out.mean<-args_simulation$out.mean[idx]
if (args_simulation$skpt_spike[idx]==0){prior.nc.skpt<-skpt_prior_default()}
if (args_simulation$skpt_spike[idx]==1){prior.nc.skpt<-skpt_prior_custom(scale=1.15)}
if (args_simulation$enth_flat[idx]==0){prior.nc.enth<-enth_prior_default()}
if (args_simulation$enth_flat[idx]==1){prior.nc.enth<-enth_prior_custom(scale=0.85)}

source("code_main.R")

Table1<-data.frame(t(outer[,,]))
Table1$idx<-idx
write.csv(Table1,file=paste0("../output/Table1/",idx,"Table1.csv"))

Table2<-data.frame(t(outer.p[,,]))
Table2$idx<-idx
write.csv(Table2,file=paste0("../output/Table2/",idx,"Table2.csv"))

Table3<-data.frame(t(outer.p.agree[,,]))
Table3$idx<-idx
write.csv(Table3,file=paste0("../output/Table3/",idx,"Table3.csv"))
#}