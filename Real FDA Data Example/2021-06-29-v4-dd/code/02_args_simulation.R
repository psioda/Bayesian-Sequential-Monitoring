#####################################
# Simulation parameters
# Evan Kwiatkowski, Feb 23, 2020
#
# The (fut|eff).mix.prob weights are
# assigned to the SKEPTICAL component
#####################################

setwd("/Users/kwiatkoe/Documents/Github/Bayesian-Sequential-Monitoring/MP_FDA_Check/code")

rm(list = ls())

# simulation1 <- expand.grid(
#   seq(.27, 0.63, length = 7), # p.IP
#   .39, # p.PC
#   2,   # freq.mntr
#   1,   # enr.shape
#   4,   # out.mean
#   0,   # fut.mix.prob
#   rep(c(seq(1, 0.25, by = -0.25)), 2), # eff.mix.prob
#   0.5,  # inf.mix.prob
#   0.05, # cred.tail
#   100,  # max.ss
#   70,   # min.ss
#   10)   # reps

# simulation2 <- expand.grid(
#   seq(.39, 0.63, length = 5), # p.IP
#   .39, # p.PC
#   2,   # freq.mntr
#   1,   # enr.shape
#   4,   # out.mean
#   0,   # fut.mix.prob
#   rep(10,152), # eff.mix.prob
#   0.5,  # inf.mix.prob
#   0.05, # cred.tail
#   100,  # max.ss
#   70,   # min.ss
#   132)   # reps
# 
# simulation <- rbind(simulation1, simulation2)

simulation <- expand.grid(
  #seq(.27, 0.63, length = 7), # p.IP
  # .39, # p.PC
  NA,
  NA,
  2,   # freq.mntr
  1,   # enr.shape
  4,   # out.mean
  0,   # fut.mix.prob
  c(seq(1, 0, by = -0.05), 10, 20, 21, 22, 23, 24, 30), # eff.mix.prob
  0.5,  # inf.mix.prob
  0.05, # cred.tail
  92,  # max.ss
  50,   # min.ss
  1)   # reps

names(simulation) <- c(
  "p.IP",
  "p.PC",
  "freq.mntr",
  "enr.shape",
  "out.mean",
  "fut.mix.prob",
  "eff.mix.prob",
  "inf.mix.prob",
  "cred.tail",
  "max.ss",
  "min.ss",
  "reps")

write.csv(x    = simulation, file = "args_simulation.csv")