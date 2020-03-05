enth_joint <- function(){
  
q.outer    <- 0.025  # y > x + delta.enth
scale      <- 1
q.inner    <- (1 - 0.8364525-q.outer)*scale   # y > x + delta.intr (major typo caught!)

source("priors/enth_marginal.R", local = TRUE)
marginal()

f1 <- function(a){   # q.outer
  enth.rd.prior <- function(x){
    exp(-(abs(x - delta.skpt)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = a[1], beta = a[2]) -
       pgnorm(q = -1, mu = delta.skpt, alpha = a[1], beta = a[2]))
  }
  integrate(enth.rd.prior, lower = delta.enth, upper = 1)$value
}
f2 <- function(a){   # q.inner
  enth.rd.prior <- function(x){
    exp(-(abs(x - delta.skpt)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))/
      (pgnorm(q = 1,  mu = delta.skpt, alpha = a[1], beta = a[2]) -
       pgnorm(q = -1, mu = delta.skpt, alpha = a[1], beta = a[2]))
  }
  integrate(enth.rd.prior, lower = delta.intr, upper = delta.enth)$value
}

start <- c(enth.rd.alpha0,enth.rd.beta0)
start
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

optim.fit <- optim(par          = start,
                   fn           = fn,
                   lower        = c(0, 0),
                   upper        = c(Inf, Inf),
                   method       = "L-BFGS-B")
optim.fit
a2 <- optim.fit$par
nlm.fit <- nlminb(start     = start,
                  objective = fn,
                  lower     = c(0, 0),
                  upper     = c(Inf, Inf))
a2 <- nlm.fit$par
a2
print(f1(a2))
print(f2(a2))
enth.rd.alpha0 <- a2[1]
enth.rd.beta0  <- a2[2]
assign("enth.rd.alpha0",enth.rd.alpha0,envir = .GlobalEnv)
assign("enth.rd.beta0",enth.rd.beta0,envir = .GlobalEnv)

source("priors/enth_conditional.R", local = TRUE)
q.outer    <- 0.1  # y > x + delta.enth
scale      <- 1
q.inner    <- (pnorm(qnorm(q.outer)/2)-q.outer)*scale   # y > x + delta.intr (major typo caught!)
mu1 <- mu + 0.2
mu2 <- mu + 0.1
enth_tail_area()

f1 <- function(a){
  enth.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
    exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
         pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
      (pgnorm(q = 1 - x, mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = 0,     mu = mu, alpha = a[1], beta  = a[2]))
  }
  enth.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
    exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
         pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
      (pgnorm(q = 1,  mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = -x, mu = mu, alpha = a[1], beta  = a[2]))
  }
  # window around mu +/- delta.enth
  c1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(enth.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  
  if (enth.prior.2(mu1,-mu1) > 1E-10){
  c3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  } else {
  c3 <- 0
  }

  c1 + c2 + c3
  }

f2 <- function(a){
  enth.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
    exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
         pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
      (pgnorm(q = 1 - x, mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = 0,     mu = mu, alpha = a[1], beta  = a[2]))
  }
  enth.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
    exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
      exp(-(abs(y - mu)/a[1])^a[2])/(2*a[1]*gamma(1/a[2])/a[2])/
      (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
         pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
      (pgnorm(q = 1,  mu = mu, alpha = a[1], beta  = a[2]) -
         pgnorm(q = -x, mu = mu, alpha = a[1], beta  = a[2]))
  }
  # window around mu +/- delta.intr
  c1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(enth.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  if (enth.prior.2(mu1,-mu1) > 1E-10){
    c3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  } else {
    c3 <- 0
  }
  d1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
  d2 <- integrate_debug(enth.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
  if (enth.prior.2(mu2,-mu2) > 1E-10){
    d3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
  } else {
    d3 <- 0
  }
  d1 + d2 + d3 - c1 - c2 - c3
  }

start <- c(enth.alpha0, enth.beta0)
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
optim.fit <- optim(par    = start,
             fn           = fn,
             lower        = c(0, 0),
             upper        = c(Inf, Inf),
             method       = "L-BFGS-B")
a2 <- optim.fit$par
a2
print(f1(a2))
print(f2(a2))
enth.alpha0 <- a2[1]
enth.beta0  <- a2[2]
assign("enth.alpha0", enth.alpha0, envir = .GlobalEnv)
assign("enth.beta0", enth.beta0, envir = .GlobalEnv)

# assemble final prior
enth.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
  exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
    exp(-(abs(y - mu)/enth.alpha0)^enth.beta0)/(2*enth.alpha0*gamma(1/enth.beta0)/enth.beta0)/
    (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
       pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
    (pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
       pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
}
enth.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
  exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
    exp(-(abs(y - mu)/enth.alpha0)^enth.beta0)/(2*enth.alpha0*gamma(1/enth.beta0)/enth.beta0)/
    (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
       pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
    (pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
       pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
}

## CHECK MARGINAL CONDITIONS
print(integrate_debug(enth.prior.1, xmin = -1,         xmax = delta.skpt, ymin = 0, ymax = function(x) 1 - x))
print(integrate_debug(enth.prior.1, xmin = delta.skpt, xmax = delta.intr, ymin = 0, ymax = function(x) 1 - x))

## CHECK CONDITIONAL CONDITIONS
c1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
c2 <- integrate_debug(enth.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
if (enth.prior.2(mu1,-mu1) > 1E-10){
  c3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
} else {
  c3 <- 0
}
d1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
d2 <- integrate_debug(enth.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
if (enth.prior.2(mu2,-mu2) > 1E-10){
  d3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
} else {
  d3 <- 0
}
print(c1 + c2 + c3)
print(d1 + d2 + d3 - c1 - c2 - c3)
}