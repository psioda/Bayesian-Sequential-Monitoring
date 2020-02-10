##################################
# Simulation parameters
# Evan Kwiatkowski, Feb 2020
##################################
rm(list = ls())

simulation <- expand.grid(
                    .39,  #seq(0.39, 0.51, by = 0.01), # p.IP
                     .39, # p.PC
                     2,   # freq.mntr
                     1,   # enr.shape
                     4,   # out.mean
                     0,   # fut.mix.prob (how much weight given to skeptical prior for futility analysis)
                     c(1,NA),   # eff.mix.prob (how much weight given to skeptical prior for efficacy analysis)
                     0.5,  # inf.mix.prob
                     0.05, # cred.tail
                     100,  # max.ss
                     1000)   # reps

names(simulation) <- c("p.IP",
                  "p.PC",
                  "freq.mntr",
                  "enr.shape",
                  "out.mean",
                  "fut.mix.prob",
                  "eff.mix.prob",
                  "inf.mix.prob",
                  "cred.tail",
                  "max.ss",
                  "reps")

write.csv(simulation, 
          file = "args_simulation.csv")

# NOTES:
# The mix.prob weights are those assigned to the SKEPTICAL component
