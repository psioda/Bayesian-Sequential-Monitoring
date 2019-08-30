rm(list = ls())

## Design parameters, usual case ##
p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
p.enth<-0.40
p.intr<-0.30
tail.skpt<-0.045  # tail probabilities for priors (low, high)
tail.enth<-0.05
sig.fut<-0.85     # significant trial result threshold
sig.eff<-0.95
cred.tail<-0.05   # credible interval is 1-cred.tail
max.ss<-76        # maximum sample size