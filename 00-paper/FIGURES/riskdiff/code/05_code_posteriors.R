##################################
### Evan Kwiatkowski, Feb 2020
##################################
# Objective: based on current data, create the following: 
# 1) scaled (un-normalized) posterior densities 
#    skpt.post.sc
#    enth.post.sc
# 2) scale factor
#    sc
# 3) scaled normalizing constant for posterior density
#    skpt.nc.sc
#    enth.nc.sc
# 4) scaled (un-normalized) marginal distributions
#    skpt.post.x.sc
#    skpt.post.y.sc
#    enth.post.x.sc
#    enth.post.y.sc
# 5) posterior mixing weights
#    fut.skpt.wt
#    fut.enth.wt
#    eff.skpt.wt
#    eff.enth.wt
#    inf.skpt.wt
#    inf.enth.wt
# Notes: Anything with "post" refers to a function

# SECTION 1: PRIORS (normalized)
skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  exp(-(abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0)/(2*skpt.rd.alpha0*gamma(1/skpt.rd.beta0)/skpt.rd.beta0)*
  exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,     mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
     pgnorm(q = -1,    mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))/
    (pgnorm(q = 1 - x, mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0) -
     pgnorm(q = 0,     mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  exp(-(abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0)/(2*skpt.rd.alpha0*gamma(1/skpt.rd.beta0)/skpt.rd.beta0)*
  exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
     pgnorm(q = -1, mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))/
    (pgnorm(q = 1,  mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0) -
     pgnorm(q = -x, mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0))
}
enth.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
  exp(-(abs(y - mu)/enth.alpha0)^enth.beta0)/(2*enth.alpha0*gamma(1/enth.beta0)/enth.beta0)/
    (pgnorm(q = 1,     mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
     pgnorm(q = -1,    mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
    (pgnorm(q = 1 - x, mu = mu,         alpha = enth.alpha0,    beta  = enth.beta0) -
     pgnorm(q = 0,     mu = mu,         alpha = enth.alpha0,    beta  = enth.beta0))
}
enth.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
    exp(-(abs(y - mu)/enth.alpha0)^enth.beta0)/(2*enth.alpha0*gamma(1/enth.beta0)/enth.beta0)/
    (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
     pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
    (pgnorm(q = 1,  mu = mu,         alpha = enth.alpha0,    beta  = enth.beta0) -
     pgnorm(q = -x, mu = mu,         alpha = enth.alpha0,    beta  = enth.beta0))
}

# SECTION 1.1 FIND MLEs (assuming nonzero cell counts in each)
PC.mle <- y1.PC/sum(y0.PC, y1.PC)
IP.mle <- y1.IP/sum(y0.IP, y1.IP)

# SECTION 2: PRIOR MIXING WEIGHTS (only if eff.mix.prob == NA)
if (is.na(eff.mix.prob)){
  if (IP.mle >= PC.mle){
  skpt.lik <- skpt.prior.1(IP.mle - PC.mle, PC.mle)
  enth.lik <- enth.prior.1(IP.mle - PC.mle, PC.mle)
  } else {
  skpt.lik <- skpt.prior.2(IP.mle - PC.mle, PC.mle)
  enth.lik <- enth.prior.2(IP.mle - PC.mle, PC.mle)
  }
  eff.mix.prob <- skpt.lik/sum(skpt.lik, enth.lik)
}

# SECTION 3: POSTERIOR DENSITIES
# log (un-normalized) posterior density
skpt.post.log.1 <- function(x, y){
  y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
  y1.PC*log(y)     + y0.PC*log(1 - y) -
  (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
  (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
  log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
      pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}
skpt.post.log.2 <- function(x, y){
  y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
  y1.PC*log(y)     + y0.PC*log(1 - y) -
  (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
  (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
  log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
      pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}
enth.post.log.1 <- function(x, y){
  y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
  y1.PC*log(y)     + y0.PC*log(1 - y) -
  (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
  (abs(y - mu)/enth.alpha0)^enth.beta0 -
  log(pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
      pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
}
enth.post.log.2 <- function(x, y){
  y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
  y1.PC*log(y)     + y0.PC*log(1 - y) -
  (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
  (abs(y - mu)/enth.alpha0)^enth.beta0 -
  log(pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
      pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
}

# scale factor: average of estimated maximum of log (un-normalized) posterior densities
if (IP.mle >= PC.mle){
  sc <- (skpt.post.log.1(IP.mle - PC.mle, PC.mle) + enth.post.log.1(IP.mle - PC.mle, PC.mle))/2
} else {
  sc <- (skpt.post.log.2(IP.mle - PC.mle, PC.mle) + enth.post.log.2(IP.mle - PC.mle, PC.mle))/2
}

# scaled (un-normalized) posterior density
skpt.post.sc.1 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
    y1.PC*log(y)     + y0.PC*log(1 - y) -
    (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
    (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
    log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
        pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
    sc
  )
}
skpt.post.sc.2 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
    y1.PC*log(y)     + y0.PC*log(1 - y) -
    (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
    (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
    log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
        pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
    sc
  )
}
enth.post.sc.1 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
    y1.PC*log(y)     + y0.PC*log(1 - y) -
    (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
    (abs(y - mu)/enth.alpha0)^enth.beta0 -
    log(pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
        pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0)) -
    sc
  )
}
enth.post.sc.2 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
    y1.PC*log(y)     + y0.PC*log(1 - y) -
    (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
    (abs(y - mu)/enth.alpha0)^enth.beta0 -
    log(pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
        pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0)) -
    sc
  )
}

# # check posterior.kernel(IP.mle - PC.mle, PC.mle)/exp(sc) = 1 (approx)
# if (IP.mle >= PC.mle){
#   print(paste0("post kernel at MLE :", 
#                (skpt.post.sc.1(IP.mle - PC.mle, PC.mle) + enth.post.sc.1(IP.mle - PC.mle, PC.mle))/2))
# } else {
#   print(paste0("post kernel at MLE :", 
#                (skpt.post.sc.2(IP.mle - PC.mle, PC.mle) + enth.post.sc.2(IP.mle - PC.mle, PC.mle))/2))
# }

# scaled normalizing constant for posterior density
skpt.nc.sc <- integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
              integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
enth.nc.sc <- integrate_debug(enth.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
              integrate_debug(enth.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

# check that posterior is normalized
#print(paste0("skpt posterior integral: ",
#      (integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#       integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
#       skpt.nc.sc))
#print(paste0("enth posterior integral: ",
#      (integrate_debug(enth.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#       integrate_debug(enth.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
#       enth.nc.sc))

# posterior mixing weights
# http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
fut.skpt.wt <- fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
eff.skpt.wt <- eff.mix.prob*skpt.nc.sc/(eff.mix.prob*skpt.nc.sc + (1 - eff.mix.prob)*enth.nc.sc)
inf.skpt.wt <- inf.mix.prob*skpt.nc.sc/(inf.mix.prob*skpt.nc.sc + (1 - inf.mix.prob)*enth.nc.sc)

# scaled (un-normalized) marginal distributions
skpt.post.x.sc.1 <- function(x, y){
  x*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
        (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
        log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
              pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
        sc
  )
}
skpt.post.x.sc.2 <- function(x, y){
  x*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
        (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
        log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
              pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
        sc
    )
}
enth.post.x.sc.1 <- function(x, y){
  x*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
        (abs(y - mu)/enth.alpha0)^enth.beta0 -
        log(pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
              pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0)) -
        sc
    )
}
enth.post.x.sc.2 <- function(x, y){
  x*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
        (abs(y - mu)/enth.alpha0)^enth.beta0 -
        log(pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
              pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0)) -
        sc
    )
}
skpt.post.y.sc.1 <- function(x, y){
  y*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
        (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
        log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
              pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
        sc
    )
}
skpt.post.y.sc.2 <- function(x, y){
  y*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
        (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
        log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
              pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
        sc
    )
}
enth.post.y.sc.1 <- function(x, y){
  y*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
        (abs(y - mu)/enth.alpha0)^enth.beta0 -
        log(pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
              pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0)) -
        sc
    )
}
enth.post.y.sc.2 <- function(x, y){
  y*
    exp(
      y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
        y1.PC*log(y)     + y0.PC*log(1 - y) -
        (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
        (abs(y - mu)/enth.alpha0)^enth.beta0 -
        log(pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
              pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0)) -
        sc
    )
}