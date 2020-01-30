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

# data 
#y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
#y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
#y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
#y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)        

y1.IP <- 25
y0.IP <- 25
y1.PC <- 25
y0.PC <- 25
fut.mix.prob <- 0
eff.mix.prob <- 1
inf.mix.prob <- 1/2

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

# scale factor: average of estimated maximum of log (un-normalized) posterior densities
sc <- 0.5*(skpt.post.log(y1.PC/sum(y1.PC,y0.PC),y1.IP/sum(y1.IP,y0.IP))+
           enth.post.log(y1.PC/sum(y1.PC,y0.PC),y1.IP/sum(y1.IP,y0.IP)))

# scaled (un-normalized) posterior density
skpt.post.sc <- function(x, y){
  exp(
    y1.PC*log(x) + y0.PC*log(1 - x) + 
      y1.IP*log(y) + y0.IP*log(1 - y) -
      (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
      (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt -
      sc
  )
}
enth.post.sc <- function(x, y){
  exp(
    y1.PC*log(x) + y0.PC*log(1 - x) + 
      y1.IP*log(y) + y0.IP*log(1 - y) -
      (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
      (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth -
      sc
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
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt -
        sc
    )
}
skpt.post.y.sc <- function(x, y){
  y*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.skpt)/sigma0.skpt)^lambda0.skpt -
        sc
    )
}
enth.post.x.sc <- function(x, y){
  x*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth -
        sc
    )
}
enth.post.y.sc <- function(x, y){
  y*
    exp(
      y1.PC*log(x) + y0.PC*log(1 - x) + 
        y1.IP*log(y) + y0.IP*log(1 - y) -
        (abs(x - mu)/sigma0.placebo)^lambda0.placebo - 
        (abs((y - x) - delta.enth)/sigma0.enth)^lambda0.enth -
        sc
    )
}
