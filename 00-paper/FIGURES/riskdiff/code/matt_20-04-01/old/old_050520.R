# un-normalized posterior density
skpt.post.1 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
            pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
  )
}
skpt.post.2 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
            pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
  )
}
enth.post.1 <- function(x, y){
  10*exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
      (abs(y - mu)/enth.alpha0)^enth.beta0 -
      log(pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
            pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
  )
}
enth.post.2 <- function(x, y){
  10*exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
      (abs(y - mu)/enth.alpha0)^enth.beta0 -
      log(pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
            pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
  )
}

# scaled normalizing constant for posterior density
skpt.nc.sc.1 <- integrate_debug(skpt.post.1, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x)
skpt.nc.sc.2 <- integrate_debug(skpt.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
enth.nc.sc.1 <- integrate_debug(enth.post.1, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x)
enth.nc.sc.2 <- integrate_debug(enth.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

skpt.nc.sc <- sum(skpt.nc.sc.1,skpt.nc.sc.2)
enth.nc.sc <- sum(enth.nc.sc.1,enth.nc.sc.2)

# posterior mixing weights
# http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
fut.skpt.wt <- fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
(1-fut.mix.prob)*enth.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
