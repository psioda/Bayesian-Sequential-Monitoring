#####################################
# Simulation parameters
# Evan Kwiatkowski, Feb 23, 2020
#
# The (fut|eff).mix.prob weights are
# assigned to the SKEPTICAL component
#####################################

# setwd("/Users/kwiatkoe/Documents/Github/Bayesian-Sequential-Monitoring/Real FDA Data Example/MP_FDA_Check/code")

rm(list = ls())

simulation <- expand.grid(
 seq(0.36,0.54,by = 0.03), # p.IP
 .39, # p.PC
 2,   # freq.mntr
 1,   # enr.shape
 4,   # out.mean
 0,   # fut.mix.prob
 c(rep(c(1, 0.5, 0), 38)), # eff.mix.prob
 # c(rep(103, 114)), # eff.mix.prob
 0.5,  # inf.mix.prob
 0.05, # cred.tail
 100,  # max.ss
 50,   # min.ss
 66)   # reps

# simulation2 <- expand.grid(
#  .63, # p.IP
#  .39, # p.PC
#  2,   # freq.mntr
#  1,   # enr.shape
#  4,   # out.mean
#  0,   # fut.mix.prob
#  c(1, 0.5, 101, 130), # eff.mix.prob
#  0.5,  # inf.mix.prob
#  0.05, # cred.tail
#  100,  # max.ss
#  50,   # min.ss
#  2)   # reps
# 
# simulation <- rbind(simulation1, simulation2)
# 
# simulation <- expand.grid(
#   #seq(.27, 0.63, length = 7), # p.IP
#   # .39, # p.PC
#   NA,
#   NA,
#   2,   # freq.mntr
#   1,   # enr.shape
#   4,   # out.mean
#   0,   # fut.mix.prob
#   # c(seq(1, 0, by = -0.05), 10, 20, 21, 22, 23, 24, 30), # eff.mix.prob
#   c(seq(1, 0, by = -0.05), 101:130), # eff.mix.prob
#   0.5,  # inf.mix.prob
#   0.05, # cred.tail
#   92,  # max.ss
#   50,   # min.ss
#   1)   # reps

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