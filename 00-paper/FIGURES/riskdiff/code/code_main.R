##################################
# Risk difference simulations
# Evan Kwiatkowski, Jan 2020
##################################

rm(list = ls())

if (.Platform$OS.type == "windows") {
  library(pracma)  # integral2
  #library(rmutil) # int2
  
  idx <- 13
}

if (.Platform$OS.type == "unix")    { 
  #library(rmutil, lib.loc = "../rpkgs/") # int2
  library(pracma, lib.loc = "../rpkgs/")  # integral2
  
  args<-commandArgs(trailingOnly = TRUE)  # sequence from batch file
  idx<-as.numeric(args[1]);
}

source("args_model.R")
source("code_integrate.R")
source("code_functions.R")
source("code_fcn_prior_placebo.R")
source("code_skpt_tail_area.R")
source("code_enth_tail_area.R")

model <- read.csv(file = "args_simulation.csv", header = TRUE, sep = ",")
p.IP <- model$p.IP[idx]
p.PC <- model$p.PC[idx]
freq.mntr <- model$freq.mntr[idx]
enr.shape <- model$enr.shape[idx]
out.mean <- model$out.mean[idx]
eff.mix.prob <- model$eff.mix.prob[idx]
fut.mix.prob <- model$fut.mix.prob[idx]
inf.mix.prob <- model$inf.mix.prob[idx]

## PRIORS ----------------------------------------------------------------
fcn_prior_placebo()
skpt_tail_area()
enth_tail_area()

## SIMULATIONS -----------------------------------------------------------
names<-c("eff.mon.initial","eff.mon.final",
                "fut.mon.initial","fut.mon.final",
                "ss.initial","ss.final",
                "post.mean.initial.IP","post.mean.final.IP",
                "post.mean.initial.PC","post.mean.final.PC",
                "cov.initial","cov.final")
inner<-array(NA, 
             dim = c(reps,length(names)), 
             dimnames = list(seq_len(reps), names))

probs.p <- c(0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1) # posterior probability
inner.p <- array(NA, 
                 dim = c(reps,2),
                 dimnames = list(seq_len(reps), c("initial.p","final.p")))

outer.p.agree <- rep(NA, 3)
names(outer.p.agree) <- c("p.agree","efficacy","conditional")

for (k in 1:reps){
  
  {print(paste0("Simulation ", k))}
  
  source("code_enrollment.R")
  
  for (h in c(seq(freq.mntr, max.ss, by = freq.mntr), max.ss)){
    
    mon.result.initial <- monitoring(index = h, 
                             fut.mix.prob = fut.mix.prob, 
                             eff.mix.prob = eff.mix.prob,
                             inf.mix.prob = inf.mix.prob)
    futility <- mon.result.initial$fut.prob
    efficacy <- mon.result.initial$eff.prob

    if (futility > sig.fut | efficacy > sig.eff){
      break
    }
  }
  
  ## INITIAL--
  n.initial <- h
  inner[k, "fut.mon.initial"] <- (futility > sig.fut)
  inner[k, "eff.mon.initial"] <- (efficacy > sig.eff)
  
  pm.cp.result.initial <- pm_cp(index = n.initial, inf.mix.prob = inf.mix.prob)
  inner[k, "post.mean.initial.PC"] <- pm.cp.result.initial$pm.mean.x
  inner[k, "post.mean.initial.IP"] <- pm.cp.result.initial$pm.mean.y
  inner[k, "cov.initial"] <- pm.cp.result.initial$coverage
  
  inner[k, "ss.initial"] <- n.initial  
  inner.p[k, "initial.p"] <- efficacy
  
  ## FINAL--
  cutoff.time <- outcome.times.all[n.initial]
  n.final <- sum(enr.times.all <= cutoff.time)
  
  mon.result.final <- monitoring(index = n.final, 
                           fut.mix.prob = fut.mix.prob, 
                           eff.mix.prob = eff.mix.prob,
                           inf.mix.prob = inf.mix.prob)
  futility.final <- mon.result.final$fut.prob
  efficacy.final <- mon.result.final$eff.prob
  inner[k,"fut.mon.final"] <- (futility.final > sig.fut)
  inner[k,"eff.mon.final"] <- (efficacy.final > sig.eff)
  
  pm.cp.result.final <- pm_cp(index = n.final, inf.mix.prob = inf.mix.prob)
  inner[k, "post.mean.final.PC"] <- pm.cp.result.final$pm.mean.x
  inner[k, "post.mean.final.IP"] <- pm.cp.result.final$pm.mean.y
  inner[k, "cov.final"] <- pm.cp.result.final$coverage
  
  inner[k,"ss.final"] <- n.final  
  inner.p[k,"final.p"] <- efficacy.final
}

outer <- apply(inner, MARGIN = 2, FUN=mean)

outer.p <- quantile(inner.p[,"final.p"][inner.p[,"initial.p"] > sig.eff & inner.p[,"final.p"] < sig.eff], probs = probs.p)

outer.p.agree["p.agree"] <- sum((inner.p[,"initial.p"] > sig.eff) == (inner.p[,"final.p"] > sig.eff)) / reps
outer.p.agree["efficacy"] <- sum((inner.p[,"initial.p"] > sig.eff) & (inner.p[,"final.p"] > sig.eff)) / reps
outer.p.agree["conditional"] <- sum((inner.p[,"initial.p"] > sig.eff) & (inner.p[,"final.p"] > sig.eff)) / sum(inner.p[,"initial.p"] > sig.eff)

Table1 <- data.frame(t(outer))
Table1$idx <- idx
write.csv(Table1, file = paste0("../output/Table1/", idx, "Table1.csv"))

Table2<-data.frame(t(outer.p))
Table2$idx <- idx
write.csv(Table2, file = paste0("../output/Table2/", idx, "Table2.csv"))

Table3<-data.frame(t(outer.p.agree))
Table3$idx <- idx
write.csv(Table3, file = paste0("../output/Table3/", idx, "Table3.csv"))