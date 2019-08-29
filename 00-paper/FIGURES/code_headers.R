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

## Design parameters, spike/slab case ##

a.s.1<-1.59
b.s.1<-6.36
a.s.2<-20
b.s.2<-80
w.s.1<-0.19
a.e.1<-3.56
b.e.1<-5.34
a.e.2<-12
b.e.2<-18
w.e.1<-0.39
y1<-5
y0<-7