##################################
### Risk difference simulations
### Evan Kwiatkowski, Feb 2020
###
### # If changes made to functions then re-run args_model.R
##################################

# setwd("/Users/kwiatkoe/Documents/Github/Bayesian-Sequential-Monitoring/Real FDA Data Example/MP_FDA_Check/code")

# for (idx in 1:8){ # check here

if (Sys.getenv("USER") %in% c("kwiatkoe", "ek50"))  {
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
set.seed(as.integer(idx*28))  #  05-19-2020

# Simulation information
simulation <- read.csv(file = "args_simulation.csv", header = TRUE, sep = ",")
for(i in 1:ncol(simulation)){
  assign(names(simulation)[i], simulation[idx, names(simulation)[i]])
}

# Simulations ---
vars  <- c("y1.IP", "y0.IP", "y1.PC", "y0.PC", "eff.prob", "fut.prob", "eff.mix.prob", "box.skpt", "box.enth", "box.ni")
suffix <- unique(c(seq(min.ss, max.ss, by = freq.mntr), max.ss))

final.inf.vars <- c("pm.skpt","pm.enth","pm.ni","pm.rd","pm.50.50","pm.skpt.enth","pm.skpt.enth.ni","pm.skpt.enth.ni.tilde","box.skpt","box.enth","box.ni","y1.IP","y1.PC","y0.IP","y0.PC")
final.inf.suffix <- "f"

names <- c(
  # paste(expand.grid(vars,suffix)$Var1, expand.grid(vars,suffix)$Var2, sep = "."), # 06-10-21 added
  paste(vars, "initial", sep = "."), 
  paste(vars, "final",   sep = "."),
  paste(final.inf.vars, final.inf.suffix, sep = "."))

inner <- array(NA, dim = c(reps, length(names)), dimnames = list(seq_len(reps), names))

for (i in 1:reps){
  
  {print(paste0("IDX ", idx, ", Simulation ", i, ", eff_mix_prob ", eff.mix.prob))}
  
  source("04_code_enrollment.R")
  
  # December 2021 - add this line
  n.initial <- max.ss
  # December 2021 - comment below
  # for (j in unique(c(seq(min.ss, max.ss, by = freq.mntr), max.ss))){
  #   #for(j in max.ss){
  #   #{print(paste0("Inner loop: ", j))}
  #   n.initial          <- j     
  #   mon.result.initial <- monitoring(index = j)
  #   futility           <- mon.result.initial$fut.prob
  #   efficacy           <- mon.result.initial$eff.prob
  #   
  #   # # 06-10-21 added
  #   # for(k in 1:length(vars)){ inner[i, paste(vars[k], j, sep = ".")] <- as.numeric(mon.result.initial[vars[k]]) }
  #   
  #   if (futility > sig.fut | efficacy > sig.eff){
  #     break
  #   }
  # }
  
  # December 2021 - comment below
  # # Initial
  # for(k in 1:length(vars)){ inner[i, paste(vars[k], "initial", sep = ".")] <- as.numeric(mon.result.initial[vars[k]]) }
  
  # Final 
  if (is.na(p.IP)){
    n.final                          <- dat[dat$targOutNum == n.initial, "nFin"]
  } else {
    cutoff.time                      <- outcome.times.all[n.initial]
    n.final                          <- sum(enr.times.all <= cutoff.time)
  }
  mon.result.final                 <- monitoring(index = n.final) # calls prior_data_conflict(), takes time
  for(k in 1:length(vars)){ inner[i, paste(vars[k], "final", sep = ".")] <- as.numeric(mon.result.final[vars[k]]) }
  
  # inf.result.final                   <- inference(index = n.final)
  # for(k in 1:length(final.inf.vars)){ inner[i, paste(final.inf.vars[k], "f", sep = ".")] <- as.numeric(inf.result.final[final.inf.vars[k]]) }
  # 
}

Table0     <- data.frame(t(inner))
Table0$idx <- idx
write.csv(Table0, file = paste0("../output/Table0/", idx, "Table0.csv"))
# }
