p.skpt     <- 0.40     # response rate for skeptic, enthusiast, futility
p.enth     <- 0.67
p.intr     <- (p.skpt+p.enth)/2
tail.skpt  <- 0.025  # tail probabilities for priors
tail.enth  <- 0.025
cred.tail  <- 0.05
sig.fut    <- 0.975
sig.eff    <- 0.975
epsilon    <- 0 # used to stop numerical error from integration
max.ss     <- 60
reps       <- 2.5E5

mu0.skpt <- p.skpt
mu0.enth <- p.enth