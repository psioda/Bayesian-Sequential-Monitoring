# normalized prior density
skpt.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
  exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
    exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
    (pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = 0,   mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
  exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
    exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
    (pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}

# tail probabilities
c1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu + delta.enth, ymax = function(x) 1 - x)
c2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu + delta.enth, ymax = 1)
result1[i,j] <- c1 + c2
d1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu + delta.intr, ymax = mu + delta.enth)
d2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu + delta.intr, ymax = mu + delta.enth)
result2[i,j] <- d1 + d2
}

integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

evan <- function(x,y) {1}
integrate_debug(evan, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
integrate_debug(evan, xmin = function(y) -y, xmax = 1, ymin = function(x) -x, ymax = function(x) 1 - x)


integrate_debug(skpt.prior.2, xmin = function(y) 1 - y, xmax = function(y) -y, ymin = 0, ymax = 1) 



# normalized prior density
skpt.prior.1 <- function(y, x){ # for x > 0 (gamma > 0)
  exp(-(abs(y - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
    exp(-(abs(x - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
    (pgnorm(q = 1 - y, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = 0,   mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}

integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
  integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

mu1<-0
integrate_debug(skpt.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x) +
integrate_debug(skpt.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1) +
integrate_debug(skpt.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
