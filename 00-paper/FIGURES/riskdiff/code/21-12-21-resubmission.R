rm(list = ls())
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-T1E-12-28-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-power-12-28-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-T1E-12-28-21.csv")
dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-power-12-28-21.csv")

summary(dat$eff.prob.final)
summary(dat$min.ss)
summary(dat$max.ss)
# summary(dat$y1.IP.initial + dat$y0.IP.initial + dat$y1.PC.initial + dat$y0.PC.initial)
summary(dat$y1.IP.final + dat$y0.IP.final + dat$y1.PC.final + dat$y0.PC.final)

mean(dat$eff.prob.final > 0.9)
mean(dat$eff.prob.final > 0.95)
mean(dat$eff.prob.final > 0.975)
mean(dat$eff.prob.final > 0.99)
# new.alpha <- 0.2335577
# new.alpha <- 0.377
# mean(dat$eff.prob.final >= 1 - new.alpha)
# mean(dat$eff.prob.final >= 0.9796)
# mean(dat$eff.prob.final >= 0.819)


library(fmsb)
res <- riskdifference(dat$y1.PC.final,
                      dat$y1.IP.final,
                      dat$y1.PC.final + dat$y0.PC.final,
                      dat$y1.IP.final + dat$y0.IP.final)
summary(res$p.value)
mean(res$p.value < 0.1)
mean(res$p.value < 0.05)
mean(res$p.value < 0.025)

# mean(res$p.value < new.alpha)

power.prop.test(n = 266, p1 = .39, p2 = .51, sig.level = 0.025, alternative = "two.sided") 
power.prop.test(n = 269, p1 = .39, p2 = .51, sig.level = 0.025, alternative = "one.sided") 
power.prop.test(n = 210, p1 = .39, p2 = .51, sig.level = 0.025, alternative = "one.sided") 

power.prop.test(n = 161, p1 = .39, p2 = .51, sig.level = 0.05, alternative = "one.sided") 

## CHANGE IN PAPER TOO