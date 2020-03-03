rm(list = ls())

library(pracma)
library(gnorm)

delta.enth <- 0.12
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2
mu         <- 0.39
sig.fut    <- 0.975
sig.eff    <- 0.975

source("skpt_joint.R")
skpt_joint()

source("enth_joint.R")
enth_joint()


