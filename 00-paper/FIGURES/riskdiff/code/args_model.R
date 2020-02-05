##################################
# Model parameters
# Evan Kwiatkowski, Feb 2020
##################################
rm(list = ls())

if (.Platform$OS.type == "windows") {
  library(pracma)
}

if (.Platform$OS.type == "unix")    { 
  library(pracma, lib.loc = "../rpkgs/")
}

delta.enth <- 0.12
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2
mu         <- 0.39
sig.fut    <- 0.975
sig.eff    <- 0.975

source("code_integrate.R")
source("code_functions.R") # contains nested source("code_posteriors.R", local = TRUE)
source("code_fcn_prior_placebo.R")
source("code_skpt_tail_area.R")
source("code_enth_tail_area.R")

fcn_prior_placebo()
skpt_tail_area()
enth_tail_area()

save.image(file = 'args_model.RData')