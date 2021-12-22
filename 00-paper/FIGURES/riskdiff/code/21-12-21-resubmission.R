# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-T1E-12-22-21.csv")
dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-power-12-22-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-T1E-12-22-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-power-12-22-21.csv")


summary(dat$eff.prob.initial)
summary(dat$eff.prob.final)
summary(dat$min.ss)
summary(dat$max.ss)
summary(dat$y1.IP.initial + dat$y0.IP.initial + dat$y1.PC.initial + dat$y0.PC.initial)
summary(dat$y1.IP.final + dat$y0.IP.final + dat$y1.PC.final + dat$y0.PC.final)

# mean(dat$eff.prob.final > 0.9)
# mean(dat$eff.prob.final >= 0.95)
mean(dat$eff.prob.final >= 0.9796)

mean(dat$eff.prob.final >= 0.819)
