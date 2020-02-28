rm(list = ls())

library(pracma)
library(gnorm)

delta.enth <- 0.12
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2
mu         <- 0.39
sig.fut    <- 0.975
sig.eff    <- 0.975
q.outer    <- 0.025  # y > x + delta.enth
q.inner    <- 0.20   # y > x + delta.intr (major typo caught!)

source("marginal.R")
marginal()

f1 <- function(a){   # q.outer
  placebo.prior <- function(x){
    exp(-(abs(x - delta.skpt)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = a[1], beta = a[2]) -
       pgnorm(q = -1, mu = delta.skpt, alpha = a[1], beta = a[2]))
  }
  integrate(placebo.prior, lower = delta.enth, upper = 1)$value
}
f2 <- function(a){   # q.inner
  placebo.prior <- function(x){
    exp(-(abs(x - delta.skpt)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = a[1], beta = a[2]) -
       pgnorm(q = -1, mu = delta.skpt, alpha = a[1], beta = a[2]))
  }
  integrate(placebo.prior, lower = delta.intr, upper = 1)$value
}

start <- c(placebo.alpha0,placebo.beta0)
start
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

nlm.fit <- nlminb(start     = start, 
                  objective = fn, 
                  lower     = c(0,0), 
                  upper     = c(Inf, Inf))
a2 <- nlm.fit$par
a2
f1(a2)
f2(a2)
placebo.alpha0 <- a2[1]
placebo.beta0  <- a2[2]

source("code_integrate.R")
source("conditional.R")
q.outer <- 0.5  # window around mu +/- delta.enth
q.inner <- 0.3  # window around mu +/- delta.intr
skpt_tail_area()

f1 <- function(a){
  skpt.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
    exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
         pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
      (pgnorm(q = 1 - x, mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = 0,     mu = mu, alpha = a[1], beta  = a[2]))
  }
  skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
    exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
         pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
      (pgnorm(q = 1,  mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = -x, mu = mu, alpha = a[1], beta  = a[2]))
  }
  # window around mu +/- delta.enth
  c1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu - delta.enth, ymax = mu + delta.enth)
  c2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu - delta.enth, ymax = mu + delta.enth)
  c1 + c2
  }

f2 <- function(a){
  skpt.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
    exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
         pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
      (pgnorm(q = 1 - x, mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = 0,     mu = mu, alpha = a[1], beta  = a[2]))
  }
  skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
    exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
         pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
      (pgnorm(q = 1,  mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = -x, mu = mu, alpha = a[1], beta  = a[2]))
  }
  # window around mu +/- delta.intr
  c1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu - delta.intr, ymax = mu + delta.intr)
  c2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu - delta.intr, ymax = mu + delta.intr)
  c1 + c2
  }

start <- c(skpt.alpha0, skpt.beta0)
start
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

# nlm.fit <- nlminb(start     = start, 
#                   objective = fn, 
#                   lower     = c(0,0), 
#                   upper     = c(Inf, Inf))
optim.fit <- optim(par    = start,
             fn           = fn,
             lower        = c(0, 0),
             upper        = c(Inf, Inf),
             method       = "L-BFGS-B")
a2 <- optim.fit$par
a2
f1(a2)
f2(a2)
skpt.alpha0 <- a2[1]
skpt.beta0  <- a2[2]

# assemble final prior
skpt.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
  exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
    exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
    (pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
  exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
    exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
    (pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}

## CHECK ORIGINAL 4 CONDITIONS
integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu - delta.intr, ymax = mu + delta.intr)+
integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu - delta.intr, ymax = mu + delta.intr)

integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu - delta.enth, ymax = mu + delta.enth)+
integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu - delta.enth, ymax = mu + delta.enth)

integrate_debug(skpt.prior.1, xmin = delta.intr, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
integrate_debug(skpt.prior.1, xmin = delta.enth, xmax = 1, ymin = 0, ymax = function(x) 1 - x)