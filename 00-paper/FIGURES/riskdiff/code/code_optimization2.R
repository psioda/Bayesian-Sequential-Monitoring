rm(list = ls())

library(pracma)

delta.enth <- 0.12
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2
mu         <- 0.39
sig.fut    <- 0.975
sig.eff    <- 0.975
q.outer    <- 0.025  # y > x + delta.enth
q.inner    <- 0.20  # y > x + delta.intr ( major typo caught! )

source("code_integrate.R")
source("code_functions.R") # contains nested source("code_posteriors.R", local = TRUE)
source("code_fcn_prior_placebo.R")
source("code_skpt_tail_area.R")

fcn_prior_placebo()
skpt_tail_area()

f1 <- function(a){
  skpt.prior <- function(x, y){ exp(-(abs(x - mu)/placebo.alpha0)^placebo.beta0 - (abs((y - x) - delta.skpt)/a[1])^a[2])}
  skpt.prior.nc <- integrate_debug(fun = skpt.prior, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  integrate_debug(fun = skpt.prior,
                  xmin = 0,
                  xmax = 1 - delta.enth,
                  ymin = function(x) x + delta.enth,
                  ymax = 1)/skpt.prior.nc
}

f2 <- function(a){
  skpt.prior <- function(x, y){ exp(-(abs(x - mu)/placebo.alpha0)^placebo.beta0 - (abs((y - x) - delta.skpt)/a[1])^a[2])}
  skpt.prior.nc <- integrate_debug(fun = skpt.prior, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  integrate_debug(fun = skpt.prior,
                  xmin = 0,
                  xmax = 1 - delta.intr,
                  ymin = function(x) x + delta.intr,
                  ymax = 1)/skpt.prior.nc
}

# start <- c(skpt.alpha0, skpt.beta0)
start <- c(0.3, 2)
start
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

nlm(fn, start)
# optim(start, fn)
a1 <- nlm(fn, start)$estimate

f1(a1)
f2(a1)

a2 <- nlminb(start     = start, 
             objective = fn, 
             lower     = c(0,0), 
             upper     = c(Inf, Inf))$par

a1
a2
