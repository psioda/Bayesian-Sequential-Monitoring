rm(list = ls())

setwd("P:/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code")

library(pracma)
library(gnorm)

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

y1.IP <- 15
y0.IP <- 5

y1.PC <- 4
y0.PC <- 16
