Table0  <- read.csv("../../output/Table0_merged.csv")

par(mfrow=c(3,2))
hist(1-Table0[(Table0$eff.mix.prob.x==10 & Table0$p.IP==0.39),"eff.mix.prob.y"], breaks = 10, freq = F, main = "", xlab = "Expected Risk Diff = 0")
hist(1-Table0[(Table0$eff.mix.prob.x==10 & Table0$p.IP==0.45),"eff.mix.prob.y"], breaks = 10, freq = F, main = "", xlab = "Expected Risk Diff = 0.06")
hist(1-Table0[(Table0$eff.mix.prob.x==10 & Table0$p.IP==0.51),"eff.mix.prob.y"], breaks = 10, freq = F, main = "", xlab = "Expected Risk Diff = 0.12")
hist(1-Table0[(Table0$eff.mix.prob.x==10 & Table0$p.IP==0.57),"eff.mix.prob.y"], breaks = 10, freq = F, main = "", xlab = "Expected Risk Diff = 0.18")
hist(1-Table0[(Table0$eff.mix.prob.x==10 & Table0$p.IP==0.63),"eff.mix.prob.y"], breaks = 10, freq = F, main = "", xlab = "Expected Risk Diff = 0.24")



table(Table0$p.IP)
