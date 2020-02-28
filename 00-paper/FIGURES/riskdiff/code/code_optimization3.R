rm(list = ls())

library(pracma)
library(gnorm)

delta.enth <- 0.12
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2
mu         <- 0.39
sig.fut    <- 0.975
sig.eff    <- 0.975

source("code_fcn_prior_placebo.R")
q.outer <- 0.5  # window around mu +/- delta.enth
q.inner <- 0.3  # window around mu +/- delta.intr
fcn_prior_placebo()
placebo.alpha0
placebo.beta0

f1 <- function(a){
  placebo.prior <- function(x){
    exp(-(abs(x-mu)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))/
      (pgnorm(q = 1, mu = mu, alpha = a[1], beta = a[2]) -
       pgnorm(q = 0, mu = mu, alpha = a[1], beta = a[2]))
  }
  integrate(placebo.prior, lower = mu - delta.enth, upper = mu + delta.enth)$value
}
f2 <- function(a){
  placebo.prior <- function(x){
    exp(-(abs(x-mu)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))/
      (pgnorm(q = 1, mu = mu, alpha = a[1], beta = a[2]) -
       pgnorm(q = 0, mu = mu, alpha = a[1], beta = a[2]))
  }
  integrate(placebo.prior, lower = mu - delta.intr, upper = mu + delta.intr)$value
}

start <- c(placebo.alpha0, placebo.beta0)
start
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

nlm.fit <- nlminb(start     = start, 
                  objective = fn, 
                  lower     = c(1E-10,1E-10), 
                  upper     = c(Inf, Inf))
a2 <- nlm.fit$par
a2
f1(a2)
f2(a2)
placebo.alpha0 <- a2[1]
placebo.beta0  <- a2[2]

source("code_integrate.R")
source("code_skpt_tail_area.R")
q.outer    <- 0.025  # y > x + delta.enth
q.inner    <- 0.20  # y > x + delta.intr ( major typo caught! )
skpt_tail_area()
skpt.alpha0
skpt.beta0

f1 <- function(a){
  skpt.prior <- function(x, y){
    exp(
      - (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
        (abs((y - x) - delta.skpt)/a[1])^a[2])/
      (pgnorm(q = 1, mu = x + delta.skpt, alpha = a[1], beta = a[2]) -
         pgnorm(q = 0, mu = x + delta.skpt, alpha = a[1], beta = a[2]))/
      (pgnorm(q = 1, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0) -
         pgnorm(q = 0, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0))*
      placebo.beta0/(2*placebo.alpha0*gamma(1/placebo.beta0))*
      a[2]/(2*a[1]*gamma(1/a[2]))
  }
  integrate_debug(fun = skpt.prior,
                  xmin = 0,
                  xmax = 1 - delta.enth,
                  ymin = function(x) x + delta.enth,
                  ymax = 1)
}

f2 <- function(a){
  skpt.prior <- function(x, y){
    exp(
      - (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
        (abs((y - x) - delta.skpt)/a[1])^a[2])/
      (pgnorm(q = 1, mu = x + delta.skpt, alpha = a[1], beta = a[2]) -
         pgnorm(q = 0, mu = x + delta.skpt, alpha = a[1], beta = a[2]))/
      (pgnorm(q = 1, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0) -
         pgnorm(q = 0, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0))*
      placebo.beta0/(2*placebo.alpha0*gamma(1/placebo.beta0))*
      a[2]/(2*a[1]*gamma(1/a[2]))
  }
  integrate_debug(fun = skpt.prior,
                  xmin = 0,
                  xmax = 1 - delta.intr,
                  ymin = function(x) x + delta.intr,
                  ymax = 1)
}

start <- c(skpt.alpha0, skpt.beta0)
start
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

nlm.fit <- nlminb(start     = start, 
                  objective = fn, 
                  lower     = c(1E-10,1E-10), 
                  upper     = c(Inf, Inf))
a2 <- nlm.fit$par
a2
f1(a2)
f2(a2)
skpt.alpha0 <- a2[1]
skpt.beta0  <- a2[2]


## check marginal distribution for theta_0
skpt.prior <- function(x, y){
  exp(
    - (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
      (abs(y - (x + delta.skpt))/skpt.alpha0)^skpt.beta0)*
    placebo.beta0/(2*placebo.alpha0*gamma(1/placebo.beta0))*
    skpt.beta0/(2*skpt.alpha0*gamma(1/skpt.beta0))/
    (pgnorm(q = 1, mu = x + delta.skpt, alpha = skpt.alpha0, beta = skpt.beta0) -
       pgnorm(q = 0, mu = x + delta.skpt, alpha = skpt.alpha0, beta = skpt.beta0))/
    (pgnorm(q = 1, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0) -
       pgnorm(q = 0, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0))
}


## CHECK ORIGINAL 4 CONDITIONS
integrate_debug(skpt.prior, mu - delta.enth, mu + delta.enth, 0, 1)
integrate_debug(skpt.prior, mu - delta.intr, mu + delta.intr, 0, 1)
integrate_debug(fun = skpt.prior,
                xmin = 0,
                xmax = 1 - delta.enth,
                ymin = function(x) x + delta.enth,
                ymax = 1)
integrate_debug(fun = skpt.prior,
                xmin = 0,
                xmax = 1 - delta.intr,
                ymin = function(x) x + delta.intr,
                ymax = 1)
