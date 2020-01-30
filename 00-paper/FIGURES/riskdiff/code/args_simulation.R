model <- expand.grid(seq(0.39, 0.51, by = 0.01), # p.IP
                     .39, # p.PC
                     2,   # freq.mntr
                     1,   # enr.shape
                     4,   # out.mean
                     0,   # fut.mix.prob (all enthusiastic)
                     1,   # eff.mix.prob (all skeptical)
                     0.5) # inf.mix.prob
names(model) <- c("p.IP",
                  "p.PC",
                  "freq.mntr",
                  "enr.shape",
                  "out.mean",
                  "fut.mix.prob",
                  "eff.mix.prob",
                  "inf.mix.prob")
write.csv(model, file = "args_simulation.csv")
