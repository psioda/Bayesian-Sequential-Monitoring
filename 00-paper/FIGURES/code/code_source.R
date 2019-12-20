rm(list = ls())

for (idx in 1:28){
if (.Platform$OS.type == "unix"){
  args<-commandArgs(trailingOnly = TRUE)
  idx<-as.numeric(args[1])}

if (.Platform$OS.type == "windows") {
  root<-"D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/code"
  setwd(root)}

source("functions.R")
source("model_args.R")
simulation_args<-read.csv(file="simulation_args.csv",header=TRUE,sep=",")

# using idx to identify model
p.range<-simulation_args$p.range[idx]
freq.mntr<-simulation_args$freq.mntr[idx]
enr.shape<-simulation_args$enr.shape[idx]
out.mean<-simulation_args$out.mean[idx]
if (simulation_args$prior_indicator_skpt[idx]==0){prior.nc.skpt<-skpt_prior_default()}
if (simulation_args$prior_indicator_skpt[idx]==1){prior.nc.skpt<-skpt_prior_custom(scale=1.15)}
if (simulation_args$prior_indicator_enth[idx]==0){prior.nc.enth<-enth_prior_default()}
if (simulation_args$prior_indicator_enth[idx]==1){prior.nc.enth<-enth_prior_custom(scale=0.85)}

source("code_main.R")

Table1<-data.frame(t(outer[,,]))
Table1$p.range<-p.range
write.csv(Table1,file=paste0("../output/Table1/",idx,"Table1.csv"))

Table2<-data.frame(t(outer.p[,,]))
Table2$p.range<-p.range
write.csv(Table2,file=paste0("../output/Table2/",idx,"Table2.csv"))

Table3<-data.frame(t(outer.p.agree[,,]))
Table3$p.range<-p.range
write.csv(Table3,file=paste0("../output/Table3/",idx,"Table3.csv"))
}