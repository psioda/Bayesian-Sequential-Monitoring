##################################
### Risk difference simulations
### Evan Kwiatkowski, Feb 2020
###
### # If changes made to functions then re-run args_model.R
###
##################################

rm(list = ls())

if (.Platform$OS.type == "windows") {
  library(pracma)
  idx <- 25
 }

if (.Platform$OS.type == "unix")    { 
  library(pracma, lib.loc = "../rpkgs/")
  args<-commandArgs(trailingOnly = TRUE)  # sequence from batch file
  idx<-as.numeric(args[1]);
}

# Model information, including all functions used. 
# The only additional source file to be called is "code_enrollment.R"
load(file = 'args_model.RData') # loads all model information include prior parameters

# Simulation information
simulation <- read.csv(file = "args_simulation.csv", header = TRUE, sep = ",")
for(i in 1:ncol(simulation)){
  assign(names(simulation)[i], simulation[idx, names(simulation)[i]])
}

# Simulations ---
names <- c("eff.mon.initial",
         "eff.mon.final",
         "fut.mon.initial",
         "fut.mon.final",
         "ss.initial",
         "ss.final",
         "mle.initial.IP",
         "mle.final.IP",
         "mle.initial.PC",
         "mle.final.PC",
         "post.mean.initial.IP",
         "post.mean.final.IP",
         "post.mean.initial.PC",
         "post.mean.final.PC",
         "cov.initial",
         "cov.final")

inner <- array(NA, 
             dim = c(reps, length(names)), 
             dimnames = list(seq_len(reps), names))

probs.p <- c(0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1) # posterior probability
inner.p <- array(NA, 
                 dim = c(reps,2),
                 dimnames = list(seq_len(reps), c("initial.p","final.p")))

outer.p.agree <- rep(NA, 3)
names(outer.p.agree) <- c("p.agree","efficacy","conditional")

for (i in 1:reps){
  
  {print(paste0("Simulation ", i))}
  
  source("code_enrollment.R")
  
  for (j in c(seq(freq.mntr, max.ss, by = freq.mntr), max.ss)){
    
    mon.result.initial <- monitoring(index = j)
    futility <- mon.result.initial$fut.prob
    efficacy <- mon.result.initial$eff.prob

    if (futility > sig.fut | efficacy > sig.eff){
      break
    }
  }
  
  # Initial ---
  n.initial <- j
  inner[i, "fut.mon.initial"] <- (futility > sig.fut)
  inner[i, "eff.mon.initial"] <- (efficacy > sig.eff)
  
  pm.cp.result.initial <- pm_cp(index = n.initial)
  inner[i, "post.mean.initial.PC"] <- pm.cp.result.initial$pm.mean.x
  inner[i, "post.mean.initial.IP"] <- pm.cp.result.initial$pm.mean.y
  inner[i, "mle.initial.PC"] <- pm.cp.result.initial$mle.PC
  inner[i, "mle.initial.IP"] <- pm.cp.result.initial$mle.IP
  inner[i, "cov.initial"] <- pm.cp.result.initial$coverage
  
  inner[i, "ss.initial"] <- n.initial  
  inner.p[i, "initial.p"] <- efficacy
  
  # Final ---
  cutoff.time <- outcome.times.all[n.initial]
  n.final <- sum(enr.times.all <= cutoff.time)
  
  mon.result.final <- monitoring(index = n.final)
  futility.final <- mon.result.final$fut.prob
  efficacy.final <- mon.result.final$eff.prob
  inner[i, "fut.mon.final"] <- (futility.final > sig.fut)
  inner[i, "eff.mon.final"] <- (efficacy.final > sig.eff)
  
  ####### Valentine's Day Debug  ####### 
  #source("valentine.R")
  #valentine(index = n.final)
  ######################################
  
  pm.cp.result.final <- pm_cp(index = n.final)
  inner[i, "post.mean.final.PC"] <- pm.cp.result.final$pm.mean.x
  inner[i, "post.mean.final.IP"] <- pm.cp.result.final$pm.mean.y
  inner[i, "cov.final"] <- pm.cp.result.final$coverage
  inner[i, "mle.final.PC"] <- pm.cp.result.final$mle.PC
  inner[i, "mle.final.IP"] <- pm.cp.result.final$mle.IP
  
  inner[i,"ss.final"] <- n.final  
  inner.p[i,"final.p"] <- efficacy.final
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