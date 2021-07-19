##################################
# Model parameters
# Evan Kwiatkowski, Feb 23, 2020
##################################
rm(list = ls())

if (Sys.getenv("USER") == "kwiatkoe") {
  library(pracma)
  library(gnorm)
  setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code")
}

if (.Platform$OS.type == "unix")    { 
  library(pracma, lib.loc = "../rpkgs/")
  library(gnorm,  lib.loc = "../rpkgs/")
}

delta.enth <- 0.12
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2
mu         <- 0.39
sig.fut    <- 0.975
sig.eff    <- 0.975

source("03_code_integrate.R")
source("06_code_functions.R") # contains nested source("code_posteriors.R", local = TRUE)

source("priors/03_enth_joint.R")
enth_joint()

source("priors/06_skpt_joint.R")
skpt_joint()

## based on skpt prior -- measures upper tail and upper half
delta.ni.enth <- 0.42   ## new upper 
delta.ni.skpt <- 0.06   ## new modal value
delta.ni.intr <- 0.24   ## new halfway point

#delta.ni.enth <- (3*delta.enth-delta.skpt)/2 ## new upper 
#delta.ni.skpt <- (delta.skpt+delta.enth)/2   ## new modal value
#delta.ni.intr <- delta.enth                  ## new halfway point
source("priors/09_ni_joint.R")
ni_joint()

save.image(file = 'args_model.RData')