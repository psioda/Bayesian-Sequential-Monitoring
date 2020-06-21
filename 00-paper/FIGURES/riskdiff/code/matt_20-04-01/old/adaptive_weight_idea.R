rm(list=ls())

binom2d <- function(x,y){
  choose(n1,x)*p1^x*(1-p1)^(n1-x)*
  choose(n2,y)*p2^y*(1-p2)^(n2-y)
}

n1 <- 10
p1 <- 0.3
n2 <- 20
p2 <- 0.6

sum <- 0
for (i in 1:(n1+1)){
  for (j in 1:(n2+1)){
    sum <- sum + binom2d(i-1,j-1)
  }
}




# Evaluates prior-data conflict metric 
# Evans2006, pg897, ex2
evans <- function(n, t0, a, b){
  
  n<-10
  a<-4
  b<-8
  t0<-5
  
  prior.kernel <- function(p) p^(a-1)*(1-p)^(b-1)
  prior.nc     <- integrate(prior.kernel, 0, 1)[[1]]


  post.nc      <- NA
  for (x in 1:(n+1)){ # has to start binomial looping at 0
    post.kernel   <- function(p) choose(n,(x-1))*p^((x-1)+a-1)*(1-p)^(n-(x-1)+b-1)/prior.nc
    post.nc[x]    <- integrate(post.kernel, 0, 1)[[1]]
  }
  
  sum(post.nc[post.nc <= post.nc[t0 + 1]])
  
}

evans(10,5,4,8)

# SECTION 1: PRIORS (normalized)
skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  dgnorm(x = x,     mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0)/
 (pgnorm(q = 1,     mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
  pgnorm(q = -1,    mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))*
  dgnorm(x = y,     mu = mu,         alpha = skpt.alpha0,    beta = skpt.beta0)/
 (pgnorm(q = 1 - x, mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0) -
  pgnorm(q = 0,     mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  dgnorm(x = x,  mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0)/
 (pgnorm(q = 1,  mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
  pgnorm(q = -1, mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))*
  dgnorm(x = y,  mu = mu,         alpha = skpt.alpha0,    beta = skpt.beta0)/
 (pgnorm(q = 1,  mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0) -
  pgnorm(q = -x, mu = mu,         alpha = skpt.alpha0,    beta  = skpt.beta0))
}


# check that posterior is normalized
integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x) +
integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)


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
  eff.mix.prob <- 0.25 + 0.75*(skpt.lik/sum(skpt.lik, enth.lik)) # 3-19-20 update
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

# check that posterior is normalized
#print(paste0("skpt posterior integral: ",
#      (integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#       integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
#       skpt.nc.sc))

# posterior mixing weights
# http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
fut.skpt.wt <- fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
eff.skpt.wt <- eff.mix.prob*skpt.nc.sc/(eff.mix.prob*skpt.nc.sc + (1 - eff.mix.prob)*enth.nc.sc)
inf.skpt.wt <- inf.mix.prob*skpt.nc.sc/(inf.mix.prob*skpt.nc.sc + (1 - inf.mix.prob)*enth.nc.sc)
