##################################
### Risk difference simulations
### Evan Kwiatkowski, Feb 2020
###
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
### # If changes made to functions then re-run args_model.R
##################################

for (idx in 23:25){ # check here

  if (Sys.getenv("USER") == "kwiatkoe") {
    library(pracma)
    library(gnorm)
  } else {                                    # longleaf
    library(pracma, lib.loc = "../rpkgs/")
    library(gnorm,  lib.loc = "../rpkgs/")
    args <- commandArgs(trailingOnly = TRUE)  # sequence from batch file
    idx  <- as.numeric(args[1]);
  }

# Model information, including all functions used (The only additional source file to be called is "code_enrollment.R")
load(file = 'args_model.RData') # loads all model information include prior parameters AND SETS SEED
set.seed(idx*623202)  #  05-19-2020

# Simulation information
simulation <- read.csv(file = "args_simulation.csv", header = TRUE, sep = ",")
for(i in 1:ncol(simulation)){
  assign(names(simulation)[i], simulation[idx, names(simulation)[i]])
}

# Simulations ---
names <- c("eff.mon.initial","eff.mon.final","fut.mon.initial","fut.mon.final",
           "ss.initial","ss.final",
           "mle.initial.IP","mle.final.IP","mle.initial.PC","mle.final.PC",
           "post.mean.initial.IP","post.mean.final.IP","post.mean.initial.PC","post.mean.final.PC",
           "cov.initial","cov.final",
           "eff.mix.prob", 
           "initial.p", "final.p",
           "box.skpt.initial", "box.skpt.final",
           "box.enth.initial", "box.enth.final",
           "risk.diff.initial", "risk.diff.final")

inner <- array(NA, 
               dim = c(reps, length(names)), 
               dimnames = list(seq_len(reps), names))

for (i in 1:reps){
  
  {print(paste0("IDX ", idx, ", Simulation ", i))}
  
  source("04_code_enrollment.R")
  
  for (j in unique(c(seq(min.ss, max.ss, by = freq.mntr), max.ss))){
    #{print(paste0("Inner loop: ", j))}
    mon.result.initial <- monitoring(index = j)
    futility           <- mon.result.initial$fut.prob
    efficacy           <- mon.result.initial$eff.prob
    if (futility > sig.fut | efficacy > sig.eff){
      break
    }
  }
  
  # Initial
  n.initial                        <- j     
  inner[i, "fut.mon.initial"]      <- (futility > sig.fut)
  inner[i, "eff.mon.initial"]      <- (efficacy > sig.eff)
  pm.cp.result.initial             <- pm_cp(index = n.initial) # calls prior_data_conflict(), takes time
  inner[i, "post.mean.initial.PC"] <- pm.cp.result.initial$pm.mean.x
  inner[i, "post.mean.initial.IP"] <- pm.cp.result.initial$pm.mean.y
  inner[i, "mle.initial.PC"]       <- pm.cp.result.initial$mle.PC
  inner[i, "mle.initial.IP"]       <- pm.cp.result.initial$mle.IP
  inner[i, "cov.initial"]          <- pm.cp.result.initial$coverage
  inner[i, "ss.initial"]           <- n.initial  
  inner[i, "initial.p"]            <- efficacy
  inner[i, "risk.diff.initial"]    <- mon.result.initial$risk.diff.mle
  inner[i, "box.skpt.initial"]     <- mon.result.initial$skpt.psi
  inner[i, "box.enth.initial"]     <- mon.result.initial$enth.psi
  
  # Final ---
  cutoff.time                      <- outcome.times.all[n.initial]
  n.final                          <- sum(enr.times.all <= cutoff.time)
  mon.result.final                 <- monitoring(index = n.final) # calls prior_data_conflict(), takes time
  futility.final                   <- mon.result.final$fut.prob
  efficacy.final                   <- mon.result.final$eff.prob
  inner[i, "fut.mon.final"]        <- (futility.final > sig.fut)
  inner[i, "eff.mon.final"]        <- (efficacy.final > sig.eff)
  inner[i, "eff.mix.prob"]         <- mon.result.final$eff.mix.prob
  pm.cp.result.final               <- pm_cp(index = n.final) # calls prior_data_conflict(), takes time
  inner[i, "post.mean.final.PC"]   <- pm.cp.result.final$pm.mean.x
  inner[i, "post.mean.final.IP"]   <- pm.cp.result.final$pm.mean.y
  inner[i, "cov.final"]            <- pm.cp.result.final$coverage
  inner[i, "mle.final.PC"]         <- pm.cp.result.final$mle.PC
  inner[i, "mle.final.IP"]         <- pm.cp.result.final$mle.IP
  inner[i,"ss.final"]              <- n.final  
  inner[i,"final.p"]               <- efficacy.final
  inner[i, "risk.diff.final"]      <- mon.result.final$risk.diff.mle
  inner[i, "box.skpt.final"]       <- mon.result.final$skpt.psi
  inner[i, "box.enth.final"]       <- mon.result.final$enth.psi
}

Table0     <- data.frame(t(inner))
Table0$idx <- idx
write.csv(Table0, file = paste0("../output/Table0/", idx, "Table0.csv"))
}