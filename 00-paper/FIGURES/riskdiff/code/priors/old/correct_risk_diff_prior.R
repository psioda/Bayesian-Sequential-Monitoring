rm(list = ls())
source("code_integrate.R")

delta.skpt     <- 0
placebo.alpha0 <- 0.5
placebo.beta0  <- 1
mu             <- 0.39
a              <- c(0.5, 1)

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
c1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x)


c2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

c1 + c2


## check marginal probs
marginal <- function(x){
  exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/
    (2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
    pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))
}
integrate(marginal, -1, 1)
integrate(marginal, 0, 0.5)
integrate_debug(skpt.prior.1, xmin = 0, xmax = 0.5, ymin = 0, ymax = function(x) 1 - x)
