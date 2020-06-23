hist(inner[,"eff.mix.prob"], freq = FALSE, 
     xlab="Skeptical Prior Weight In Mixture",
     main="Skeptical Prior Weight In Mixture")

plot(inner[,"mle.final.IP"]-inner[,"mle.final.PC"], inner[,"eff.mix.prob"],
     xlab = "Risk Difference",
     ylab = "Skeptical Prior Weight In Mixture",
     main = "Skeptical Weight by Risk Difference")
plot(inner[,"ss.initial"], inner[,"eff.mix.prob"],
     xlab = "Sample Size at Enrollment Termination",
     ylab = "Skeptical Prior Weight in Mixture",
     main = "Skeptical Weight by Sample Size")
lines(lowess(inner[,"ss.initial"], inner[,"eff.mix.prob"]))

