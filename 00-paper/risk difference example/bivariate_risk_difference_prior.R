rm(list = ls())

if (.Platform$OS.type == "windows") {
  root<-"D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/risk difference example"
  setwd(root)
  
  require(rmutil)
  require(lattice)
  require(pracma)
  
  source("bivariate_functions.R")
  idx<-1
  
  # prior parameters
  sigma<-8
  alpha<-2
  delta.enth<-0.1
  delta.skpt<-0
  delta.intr=(delta.skpt+delta.enth)/2
  mu<-0.2
  
  # simulation parameters
  sig.fut<-0.85
  sig.eff<-0.95
  cred.tail<-0.05
  max.ss<-100
  reps<-10
  p.IP<-0.2
  p.PC<-0.2
  freq.mntr<-10
  enr.shape<-1
  out.mean<-4
}

if (.Platform$OS.type == "unix")    { 
  # sequence from batch file
  args<-commandArgs(trailingOnly = TRUE)
  idx<-as.numeric(args[1]);
  
  # set directories
  root<-"/nas/longleaf/home/ekwiatko/FDA/riskdiff/"
  rpkgs<-paste0(root,"rkpgs")
  output<-paste0(root,"output/")
  setwd(output)
  
  # load packages
  require(rmutil,lib.loc="/nas/longleaf/home/ekwiatko/FDA/riskdiff/rpkgs/")
  require(lattice)
  require(pracma,lib.loc="/nas/longleaf/home/ekwiatko/FDA/riskdiff/rpkgs/")
  
  # source functions
  source("/nas/longleaf/home/ekwiatko/FDA/riskdiff/code/functions.R")
  source("/nas/longleaf/home/ekwiatko/FDA/riskdiff/code/parms.R")
  
  # using idx to identify model
  model<-read.csv(file="/nas/longleaf/home/ekwiatko/FDA/riskdiff/code/model.csv",header=TRUE,sep=",")
  p.IP<-model$p.IP[idx]
  p.PC<-model$p.PC[idx]
  freq.mntr<-model$freq.mntr[idx]
  enr.shape<-model$enr.shape[idx]
  out.mean<-model$out.mean[idx]
}

## SIMULATIONS ############################################################

names<-c("eff.mon.initial","eff.mon.final",
                "fut.mon.initial","fut.mon.final",
                "ss.initial","ss.final",
                "post.mean.initial.IP","post.mean.final.IP",
                "post.mean.initial.PC","post.mean.final.PC",
                "cov.initial","cov.final")

inner<-array(NA,dim=c(reps,length(names)),dimnames=list(seq_len(reps),names))

# posterior probability
probs.p<-c(0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1)
inner.p<-array(NA, dim=c(reps,2),dimnames=list(seq_len(reps),c("initial.p","final.p")))
outer.p.agree<-rep(NA,3)
names(outer.p.agree)<-c("p.agree","efficacy","conditional")


for (k in 1:reps){
  
  if (k%%10==0){print(paste0("Simulation ",k))}
  
  group<-c(sample(c(rep("PC",4),rep("IP",20))),sample(c(rep("PC",38),rep("IP",38))))
  enr.times.all<-seq(1:100)*17
  outcome.times.all<-enr.times.all+52
  
  enr.times.PC<-enr.times.all[group=="PC"]
  outcome.times.PC<-outcome.times.all[group=="PC"]
  responses.PC<-rbinom(n=42,size=1,prob=p.PC)
  
  enr.times.IP<-enr.times.all[group=="IP"]
  outcome.times.IP<-outcome.times.all[group=="IP"]
  responses.IP<-rbinom(n=58,size=1,prob=p.IP)
  
  # outcome.times.all<-sort(c(outcome.times.IP,outcome.times.PC))
  # enr.times.all<-sort(c(enr.times.IP,enr.times.PC))
  # 
  # enr.times.IP<-cumsum(rgamma(n=max.ss,shape=enr.shape,scale=0.5))
  # outcome.times.IP<-sort(enr.times.IP+rnorm(n=max.ss,mean=out.mean,sd=0.25))
  # responses.IP<-rbinom(n=max.ss,size=1,prob=p.IP)
  # 
  # enr.times.PC<-cumsum(rgamma(n=max.ss,shape=enr.shape,scale=0.5))
  # outcome.times.PC<-sort(enr.times.PC+rnorm(n=max.ss,mean=out.mean,sd=0.25))
  # responses.PC<-rbinom(n=max.ss,size=1,prob=p.PC)
  

  
  for (h in seq(freq.mntr,(max.ss),by=freq.mntr)){
    result<-eff_fut(h)
    futility<-result[2]
    efficacy<-result[1]
    if (futility>sig.fut | efficacy>sig.eff){
      break
    }
  }
  
  ## INITIAL ##
  n.initial<-h
  result.initial<-pm_cp(n.initial)
  inner[k,"fut.mon.initial"]<-(futility>sig.fut)
  inner[k,"eff.mon.initial"]<-(efficacy>sig.eff)
  inner[k,"post.mean.initial.PC"]<-result.initial[1]
  inner[k,"post.mean.initial.IP"]<-result.initial[2]
  inner[k,"cov.initial"]<-result.initial[3]
  inner[k,"ss.initial"]<-n.initial  
  inner.p[k,"initial.p"]<-efficacy
  
  ## FINAL ##
  cutoff.time<-outcome.times.all[n.initial]
  n.final<-sum(enr.times.all<=cutoff.time)
  result<-eff_fut(n.final)
  futility.final<-result[2]
  efficacy.final<-result[1]
  result.final<-pm_cp(n.final)
  
  inner[k,"fut.mon.final"]<-(futility.final>sig.fut)
  inner[k,"eff.mon.final"]<-(efficacy.final>sig.eff)
  inner[k,"post.mean.final.PC"]<-result.final[1]
  inner[k,"post.mean.final.IP"]<-result.final[2]
  inner[k,"cov.final"]<-result.final[3]
  inner[k,"ss.final"]<-n.final  
  inner.p[k,"final.p"]<-efficacy.final
}

outer<-apply(inner,MARGIN=2,FUN=mean)

outer.p<-quantile(inner.p[,"final.p"][inner.p[,"initial.p"]>sig.eff & inner.p[,"final.p"]<sig.eff],
                        probs=probs.p)

outer.p.agree["p.agree"]<-
  sum((inner.p[,"initial.p"]>sig.eff)==(inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree["efficacy"]<-
  sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree["conditional"]<-
  sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/sum(inner.p[,"initial.p"]>sig.eff)

# output files
if (.Platform$OS.type == "windows")    { 
  write.csv(outer,file=paste0(idx,"Table1.csv"))
  write.csv(outer.p,file=paste0(idx,"Table2.csv"))
  write.csv(outer.p.agree,file=paste0(idx,"Table3.csv"))
}
if (.Platform$OS.type == "unix")    { 
write.csv(outer,file=paste0(output,"Table1/",idx,"Table1.csv"))
write.csv(outer.p,file=paste0(output,"Table2/",idx,"Table2.csv"))
write.csv(outer.p.agree,file=paste0(output,"Table3/",idx,"Table3.csv"))
}