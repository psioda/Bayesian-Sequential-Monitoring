##################################
# Model parameters
# Evan Kwiatkowski, Feb 23, 2020
##################################
rm(list = ls())

if (.Platform$OS.type == "windows") {
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

source("code_integrate.R")
source("code_functions.R") # contains nested source("code_posteriors.R", local = TRUE)

source("priors/skpt_joint.R")
skpt_joint()

source("priors/enth_joint.R")
enth_joint()

save.image(file = 'args_model.RData')
