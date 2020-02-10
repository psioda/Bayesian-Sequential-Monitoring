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

# log (un-normalized) posterior density
skpt.post.log <- function(x, y){
  y1.PC*log(x) + y0.PC*log(1 - x) + 
    y1.IP*log(y) + y0.IP*log(1 - y) -
    (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
    (abs((y - x) - delta.skpt)/skpt.alpha0)^skpt.beta0
}
enth.post.log <- function(x, y){
  y1.PC*log(x) + y0.PC*log(1 - x) + 
    y1.IP*log(y) + y0.IP*log(1 - y) -
    (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
    (abs((y - x) - delta.enth)/enth.alpha0)^enth.beta0
}

# scale factor: average of estimated maximum of log (un-normalized) posterior densities
if (min(y1.IP,y0.IP,y1.PC,y0.PC) > 0){
sc <- 0.5*(skpt.post.log(y1.PC/sum(y1.PC, y0.PC), y1.IP/sum(y1.IP, y0.IP)) + 
             enth.post.log(y1.PC/sum(y1.PC, y0.PC), y1.IP/sum(y1.IP, y0.IP)))
} else {
sc <- 0.5*(skpt.post.log(p.PC,p.IP) + 
             enth.post.log(p.PC,p.IP))
}

# scaled (un-normalized) posterior density
skpt.post.sc <- function(x, y){
  exp(
    y1.PC*log(x) + y0.PC*log(1 - x) + 
      y1.IP*log(y) + y0.IP*log(1 - y) -
      (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
      (abs((y - x) - delta.skpt)/skpt.alpha0)^skpt.beta0 -
      sc
  )
}
enth.post.sc <- function(x, y){
  exp(
    y1.PC*log(x) + y0.PC*log(1 - x) + 
      y1.IP*log(y) + y0.IP*log(1 - y) -
      (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
      (abs((y - x) - delta.enth)/enth.alpha0)^enth.beta0 -
      sc
  )
}

# scaled normalizing constant for posterior density
skpt.nc.sc <- integrate_debug(fun = skpt.post.sc,
                                xmin = 0,
                                xmax = 1,
                                ymin = 0,
                                ymax = 1)
enth.nc.sc <- integrate_debug(fun = enth.post.sc,
                                xmin = 0,
                                xmax = 1,
                                ymin = 0,
                                ymax = 1)

# # prior mixing weights (Feb 10, 2020)
if (is.na(eff.mix.prob)){
  if (min(y1.IP,y0.IP,y1.PC,y0.PC) > 0){
    PC.mle <- y1.PC/sum(y0.PC, y1.PC)
    IP.mle <- y1.IP/sum(y0.IP, y1.IP)
  }
  else {
    PC.mle <- p.PC
    IP.mle <- p.IP
    }
  skpt.lik <- skpt.post.sc(PC.mle, IP.mle)
  enth.lik <- enth.post.sc(PC.mle, IP.mle)
  eff.mix.prob <- skpt.lik/sum(skpt.lik, enth.lik)
}

# posterior mixing weights
# http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
fut.skpt.wt <- fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
eff.skpt.wt <- eff.mix.prob*skpt.nc.sc/(eff.mix.prob*skpt.nc.sc + (1 - eff.mix.prob)*enth.nc.sc)
inf.skpt.wt <- inf.mix.prob*skpt.nc.sc/(inf.mix.prob*skpt.nc.sc + (1 - inf.mix.prob)*enth.nc.sc)

# scaled (un-normalized) marginal distributions
skpt.post.x.sc <- function(x, y){
  x*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
        (abs((y - x) - delta.skpt)/skpt.alpha0)^skpt.beta0 -
        sc
    )
}
skpt.post.y.sc <- function(x, y){
  y*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
        (abs((y - x) - delta.skpt)/skpt.alpha0)^skpt.beta0 -
        sc
    )
}
enth.post.x.sc <- function(x, y){
  x*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
        (abs((y - x) - delta.enth)/enth.alpha0)^enth.beta0 -
        sc
    )
}
enth.post.y.sc <- function(x, y){
  y*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
        (abs((y - x) - delta.enth)/enth.alpha0)^enth.beta0 -
        sc
    )
}

