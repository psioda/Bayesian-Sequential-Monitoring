## Posterior Mean & Coverage Probability ------------
pm_cp <- function(index, inf.mix.prob){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)
  
  source("code_posteriors.R", local = TRUE)
  
  # posterior means
  skpt.pm.x <- (integrate_debug(fun = skpt.post.x.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1))/skpt.nc.sc
  skpt.pm.y <- (integrate_debug(fun = skpt.post.y.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1))/skpt.nc.sc
  enth.pm.x <- (integrate_debug(fun = enth.post.x.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1))/enth.nc.sc
  enth.pm.y <- (integrate_debug(fun = enth.post.y.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1))/enth.nc.sc
  
  # posterior means weighted
  pm.mean.x <- inf.skpt.wt*skpt.pm.x + (1 - inf.skpt.wt)*enth.pm.x
  pm.mean.y <- inf.skpt.wt*skpt.pm.y + (1 - inf.skpt.wt)*enth.pm.y
  
  ## coverage probability
  # create grid and find index of point representing true response probabilities
  grid.x <- p.PC
  grid.y <- p.IP
  grid <- expand.grid(seq(0,1,by=0.01),seq(0,1,by=0.01))
  grid.index <- which.min(sqrt((grid$Var1-grid.x)^2+(grid$Var2-grid.y)^2))
  
  # evaluate "normalized density" at the grid
  grid.eval <- (inf.skpt.wt/skpt.nc.sc*skpt.post.sc(grid$Var1,grid$Var2) +
    (1-inf.skpt.wt)/enth.nc.sc*enth.post.sc(grid$Var1,grid$Var2))/
    sum(inf.skpt.wt/skpt.nc.sc*skpt.post.sc(grid$Var1,grid$Var2) +
        (1-inf.skpt.wt)/enth.nc.sc*enth.post.sc(grid$Var1,grid$Var2))
  
  # grid.eval[is.nan(grid.eval)]  <-  0 # 10-29-2019, only if necessary
  
  # find the cutoff for cred tail percentile (absolute value unnecessary)
  grid.index2 <- which.min(abs(cumsum(sort(grid.eval)) - cred.tail))
  
  # is grid point in credible interval?
  coverage <- grid.eval[grid.index] >= sort(grid.eval)[grid.index2]  
  
  return(data.frame(cbind(pm.mean.x,pm.mean.y,coverage)))
}

## Monitoring ----------------------------------------------------------
monitoring <- function(index, fut.mix.prob, eff.mix.prob, inf.mix.prob){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)
  
  source("code_posteriors.R", local = TRUE)
  
  # efficacy and futility probabilities
  eff.prob.skpt <- integrate_debug(skpt.post.sc, 
                             xmin = 0, 
                             xmax = 1 - delta.skpt,
                             ymin = function(x) x + delta.skpt,
                             ymax = 1)
  eff.prob.enth <- integrate_debug(enth.post.sc, 
                             xmin = 0, 
                             xmax = 1 - delta.skpt,
                             ymin = function(x) x + delta.skpt,
                             ymax = 1)
  fut.prob.skpt <- 1 - integrate_debug(skpt.post.sc, 
                                 xmin = 0,
                                 xmax = 1 - delta.intr,
                                 ymin = function(x) x + delta.intr,
                                 ymax = 1)
  fut.prob.enth <- 1 - integrate_debug(enth.post.sc, 
                                 xmin = 0, 
                                 xmax = 1-delta.intr,
                                 ymin = function(x) x + delta.intr,
                                 ymax = 1)
  
  eff.prob <- (eff.skpt.wt/skpt.nc.sc)*eff.prob.skpt + ((1 - eff.skpt.wt)*enth.nc.sc)*eff.prob.enth
  inf.prob <- (inf.skpt.wt/skpt.nc.sc)*eff.prob.skpt + ((1 - inf.skpt.wt)*enth.nc.sc)*eff.prob.enth
  fut.prob <- (fut.skpt.wt/skpt.nc.sc)*fut.prob.skpt + ((1 - fut.skpt.wt)*enth.nc.sc)*fut.prob.enth
  
  return(data.frame(cbind(eff.prob, inf.prob, fut.prob)))
}
