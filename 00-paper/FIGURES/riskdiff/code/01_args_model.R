##################################
# Model parameters
# Evan Kwiatkowski, Feb 23, 2020
##################################
rm(list = ls())

if (Sys.getenv("USER") == "kwiatkoe") {
  library(pracma)
  library(gnorm)
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
save.image(file = 'args_model.RData')