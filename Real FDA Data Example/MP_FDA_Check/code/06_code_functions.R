##################################
### Evan Kwiatkowski, Feb 2020
###
### contains nested source("code_posteriors.R", local = TRUE)
###
##################################

# Monitoring ---
monitoring <- function(index){
  
  if (is.na(p.IP)){
    y1.IP <- dat[dat$targOutNum == index, "yObs1"]
    y0.IP <- dat[dat$targOutNum == index, "nObs1"] - dat[dat$targOutNum == index, "yObs1"]
    y1.PC <- dat[dat$targOutNum == index, "yObs0"]
    y0.PC <- dat[dat$targOutNum == index, "nObs0"] - dat[dat$targOutNum == index, "yObs0"]
  } else {
    y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
    y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
    y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
    y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)
  }
  source("05_code_posteriors.R", local = TRUE)
  
  ## efficacy probability using skeptical prior (scaled)
  #  integrate joint prior from [delta.skpt, 1], stop if greater than 1 - epsilon
  eff.prob.skpt <- integrate_debug(skpt.post.sc.1, xmin = delta.skpt, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
  
  ## efficacy probability using enthusiastic prior (scaled)
  #  integrate joint prior from [delta.skpt, 1], stop if greater than 1 - epsilon
  eff.prob.enth <- integrate_debug(enth.post.sc.1, xmin = delta.skpt, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
  
  ## futility probability using skeptical prior (scaled)
  # integrate joint prior from [delta.intr , 1], stop if less than epsilon
  # alternatively, stop if 1 - (fut.prob.skpt/skpt.nc.sc) is greater than 1 - epsilon
  # UPDATE 9/29/20: change xmin from delta.intr to delta.enth
  fut.prob.skpt <- integrate_debug(skpt.post.sc.1, xmin = delta.enth, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
  
  # futility probability using enthusiastic prior (scaled)
  # integrate joint prior from [delta.intr , 1], stop if less than epsilon
  # alternatively, stop if 1 - (fut.prob.skpt/skpt.nc.sc) is greater than 1 - epsilon
  # UPDATE 9/29/20: change xmin from delta.intr to delta.enth
  fut.prob.enth <- integrate_debug(enth.post.sc.1, xmin = delta.enth, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
  
  # posterior mixing weights
  # http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
  fut.skpt.wt <- fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
  eff.skpt.wt <- eff.mix.prob*skpt.nc.sc/(eff.mix.prob*skpt.nc.sc + (1 - eff.mix.prob)*enth.nc.sc)
  inf.skpt.wt <- inf.mix.prob*skpt.nc.sc/(inf.mix.prob*skpt.nc.sc + (1 - inf.mix.prob)*enth.nc.sc)
  
  # Recall default is eff.skpt.wt <- 1   (all skeptical prior for efficacy monitoring)
  eff.prob <- eff.skpt.wt*(eff.prob.skpt/skpt.nc.sc) + (1 - eff.skpt.wt)*(eff.prob.enth/enth.nc.sc)
  
  # Recall default is fut.skpt.wt <- 0   (all enthusiastic prior for futility monitoring)
  fut.prob <- fut.skpt.wt*(1 - fut.prob.skpt/skpt.nc.sc) + (1 - fut.skpt.wt)*(1 - fut.prob.enth/enth.nc.sc)
  
  # Recall inference prior is also making a judgement of efficacy
  # Recall default is inf.skpt.wt <- 0.5 (50:50 mixture for inference)
  inf.prob <- inf.skpt.wt*(eff.prob.skpt/skpt.nc.sc) + (1 - inf.skpt.wt)*(eff.prob.enth/enth.nc.sc)
  
  # bug caught 2/1/20, (1 - fut.prob.skpt)) must be computed before dividing by normalizing constant
  # bug caught 2/2/20, 1 - fut.prob.enth/enth.nc.sc NOT (1 - fut.prob.enth)/enth.nc.sc
  # bug caught 3/3/20, first term of fut.prob was fut.skpt.wt*((1 - fut.prob.enth/enth.nc.sc)
                                      # should be fut.skpt.wt*((1 - fut.prob.skpt/skpt.nc.sc)
  
  return(data.frame(cbind(eff.prob, inf.prob, fut.prob, eff.mix.prob, risk.diff.mle, box.skpt, box.enth, box.ni, y1.IP, y1.PC, y0.IP, y0.PC)))
}