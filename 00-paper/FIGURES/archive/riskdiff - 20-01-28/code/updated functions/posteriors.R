# Objective: based on current data, create the following: 
# 1) scaled (un-normalized) posterior densities 
#    skpt.post.sc
#    enth.post.sc
# 2) scale factor
#    skpt.sc
#    enth.sc
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

# data 
y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)        

# log (un-normalized) posterior density
skpt.post.log <- function(x, y){
  y1.PC*log(x) + y0.PC*log(1 - x) + 
    y1.IP*log(y) + y0.IP*log(1 - y) -
    (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
    (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt
}
enth.post.log <- function(x, y){
  y1.PC*log(x) + y0.PC*log(1 - x) + 
    y1.IP*log(y) + y0.IP*log(1 - y) -
    (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
    (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth
}

# scale factor: estimated maximum of log (un-normalized) posterior density
skpt.sc <- skpt.post.log(y1.PC/sum(y1.PC,y0.PC),y1.IP/sum(y1.IP,y0.IP))
enth.sc <- enth.post.log(y1.PC/sum(y1.PC,y0.PC),y1.IP/sum(y1.IP,y0.IP))

# scaled (un-normalized) posterior density
skpt.post.sc <- function(x, y){
  exp(
    y1.PC*log(x) + y0.PC*log(1 - x) + 
      y1.IP*log(y) + y0.IP*log(1 - y) -
      (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
      (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt -
      skpt.sc
  )
}
enth.post.sc <- function(x, y){
  exp(
    y1.PC*log(x) + y0.PC*log(1 - x) + 
      y1.IP*log(y) + y0.IP*log(1 - y) -
      (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
      (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth -
      enth.sc
  )
}

# scaled normalizing constant for posterior density
skpt.nc.sc <- integral2(fun = skpt.post.sc,
                                xmin = 0,
                                xmax = 1,
                                ymin = 0,
                                ymax = 1)$Q
enth.nc.sc <- integral2(fun = enth.post.sc,
                                xmin = 0,
                                xmax = 1,
                                ymin = 0,
                                ymax = 1)$Q

# adding on log scale 
# https://stackoverflow.com/questions/5802592/dealing-with-very-small-numbers-in-r
logxpy <- function(lx,ly) max(lx,ly) + log1p(exp(-abs(lx-ly)))

# posterior mixing weights
# http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
fut.skpt.wt <- exp(log(fut.mix.prob) + log(skpt.nc.sc) + skpt.sc -
                     logxpy(log(fut.mix.prob) + log(skpt.nc.sc) + skpt.sc,
                            log(1 - fut.mix.prob) + log(enth.nc.sc) + enth.sc))
fut.enth.wt <- exp(log(1 - fut.mix.prob) + log(enth.nc.sc) + enth.sc -
                     logxpy(log(fut.mix.prob) + log(skpt.nc.sc) + skpt.sc,
                            log(1 - fut.mix.prob) + log(enth.nc.sc) + enth.sc))
eff.skpt.wt <- exp(log(eff.mix.prob) + log(skpt.nc.sc)+skpt.sc -
                     logxpy(log(eff.mix.prob) + log(skpt.nc.sc) + skpt.sc,
                            log(1 - eff.mix.prob) + log(enth.nc.sc) + enth.sc))
eff.enth.wt <- exp(log(1 - eff.mix.prob) + log(enth.nc.sc) + enth.sc -
                     logxpy(log(eff.mix.prob)+log(skpt.nc.sc) + skpt.sc,
                            log(1 - eff.mix.prob)+log(enth.nc.sc) + enth.sc))
inf.skpt.wt <- exp(log(inf.mix.prob) + log(skpt.nc.sc) + skpt.sc -
                     logxpy(log(inf.mix.prob)+log(skpt.nc.sc) + skpt.sc,
                            log(1 - inf.mix.prob)+log(enth.nc.sc) + enth.sc))
inf.enth.wt <- exp(log(1 - inf.mix.prob) + log(enth.nc.sc) + enth.sc -
                     logxpy(log(inf.mix.prob) + log(skpt.nc.sc) + skpt.sc,
                            log(1 - inf.mix.prob) + log(enth.nc.sc) + enth.sc))

# scaled (un-normalized) marginal distributions
skpt.post.x.sc <- function(x, y){
  x*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt -
        skpt.sc
    )
}
skpt.post.y.sc <- function(x, y){
  y*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt -
        skpt.sc
    )
}
enth.post.x.sc <- function(x, y){
  x*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth -
        enth.sc
    )
}
enth.post.y.sc <- function(x, y){
  y*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth -
        enth.sc
    )
}