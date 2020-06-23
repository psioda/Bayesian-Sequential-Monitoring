#####################################
# Simulation parameters
# Evan Kwiatkowski, Feb 23, 2020
#
# The (fut|eff).mix.prob weights are
# assigned to the SKEPTICAL component
#####################################

rm(list = ls())

# simulation <- expand.grid(
#   seq(.39, 0.63, length = 2), # p.IP
#   .39, # p.PC
#   2,   # freq.mntr
#   1,   # enr.shape
#   4,   # out.mean
#   0,   # fut.mix.prob
#   seq(0.25,1,by=0.25),  # eff.mix.prob
#   0.5,  # inf.mix.prob
#   0.05, # cred.tail
#   100,  # max.ss
#   70,   # min.ss
#   100)   # reps
# 
# simulation <- expand.grid(
#                    seq(.39, 0.63, length = 25), # p.IP
#                     .39, # p.PC
#                     2,   # freq.mntr
#                     1,   # enr.shape
#                     4,   # out.mean
#                     0,   # fut.mix.prob
#                     c(seq(1, 0.25, by = -0.25), NA), # eff.mix.prob
#                     0.5,  # inf.mix.prob
#                     0.05, # cred.tail
#                     100,  # max.ss
#                     70,   # min.ss
#                     50)   # reps

# simulation1 <- expand.grid(
#                    seq(.39, 0.63, length = 5), # p.IP
#                     .39, # p.PC
#                     10,   # freq.mntr
#                     1,   # enr.shape
#                     4,   # out.mean
#                     0,   # fut.mix.prob
#                     c(seq(1, 0.25, by = -0.25), NA), # eff.mix.prob
#                     0.5,  # inf.mix.prob
#                     0.05, # cred.tail
#                     100,  # max.ss
#                     70,   # min.ss
#                     1000)   # reps

simulation2 <- expand.grid(
                   seq(.39, 0.63, length = 9), # p.IP
                    rep(.39, 20), # p.PC
                    2,   # freq.mntr
                    1,   # enr.shape
                    4,   # out.mean
                    0,   # fut.mix.prob
                    c(seq(1,0.25,by=-0.25)),  # eff.mix.prob (10 means prior data conflcit)
                    0.5,  # inf.mix.prob
                    0.05, # cred.tail
                    100,  # max.ss
                    70,   # min.ss
                    1)   # reps

#simulation <- rbind(simulation1, simulation2)
simulation  <- simulation2

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