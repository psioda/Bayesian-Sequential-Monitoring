rm(list = ls())

if (.Platform$OS.type == "windows") {

  library(rmutil)
  library(lattice)
  library(pracma)
  idx <- 52 # 1 14 27 40 53
}

if (.Platform$OS.type == "unix")    { 
  args<-commandArgs(trailingOnly = TRUE)   # sequence from batch file
  idx<-as.numeric(args[1]);
  
  library(rmutil, lib.loc = "../rpkgs/")
  library(lattice)
  library(pracma, lib.loc = "../rpkgs/")
}

source("code_functions.R")
source("args_model.R")

model <- read.csv(file = "args_simulation.csv", header = TRUE, sep = ",")
p.IP <- model$p.IP[idx]
p.PC <- model$p.PC[idx]
freq.mntr <- model$freq.mntr[idx]
enr.shape <- model$enr.shape[idx]
out.mean <- model$out.mean[idx]
mix.prob <- model$mix.prob[idx]

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
  
  if (k%%10 == 0){print(paste0("Simulation ", k))}
  
  group <- c(sample(c(rep("PC",4), rep("IP",20))),
             sample(c(rep("PC",38), rep("IP",38))))
  enr.times.all <- seq(1:100) * 17
  outcome.times.all <- enr.times.all + 52
  
  enr.times.PC <- enr.times.all[group == "PC"]
  outcome.times.PC <- outcome.times.all[group == "PC"]
  responses.PC <- rbinom(n = 42, size = 1, prob = p.PC)
  
  enr.times.IP <- enr.times.all[group == "IP"]
  outcome.times.IP <- outcome.times.all[group == "IP"]
  responses.IP <- rbinom(n= 58, size = 1, prob = p.IP)
  
  for (h in seq(freq.mntr, max.ss, by = freq.mntr)){ # need to add max.ss if freq.mntr not a factor
    
    result <- eff_fut(h)
    futility <- result[2]
    efficacy <- inference.skpt(h)

    if (futility > sig.fut | efficacy > sig.eff){
      break
    }
  }
  
  ## INITIAL--
  n.initial <- h
  result.initial <- pm_cp(n.initial)
  inner[k, "fut.mon.initial"] <- (futility > sig.fut)
  inner[k, "eff.mon.initial"] <- (efficacy > sig.eff)
  inner[k, "post.mean.initial.PC"] <- result.initial[1]
  inner[k, "post.mean.initial.IP"] <- result.initial[2]
  inner[k, "cov.initial"] <- result.initial[3]
  inner[k, "ss.initial"] <- n.initial  
  inner.p[k, "initial.p"] <- efficacy
  
  ## FINAL--
  cutoff.time <- outcome.times.all[n.initial]
  n.final <- sum(enr.times.all <= cutoff.time)
  
  result <- eff_fut(n.final)
  futility.final <- result[2]
  efficacy.final <- inference.skpt(n.final)
  result.final <- pm_cp(n.final)
  
  inner[k,"fut.mon.final"] <- (futility.final > sig.fut)
  inner[k,"eff.mon.final"] <- (efficacy.final > sig.eff)
  inner[k,"post.mean.final.PC"] <- result.final[1]
  inner[k,"post.mean.final.IP"] <- result.final[2]
  inner[k,"cov.final"] <- result.final[3]
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