## Posterior Mean & Coverage Probability ------------
pm_cp <- function(index, inf.mix.prob){
  
  source(monitoring.R)
  
  skpt.pm.x <- (integral2(fun = skpt.post.x.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q)/skpt.nc.sc
  skpt.pm.y <- (integral2(fun = skpt.post.y.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q)/skpt.nc.sc
  enth.pm.x <- (integral2(fun = enth.post.x.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q)/enth.nc.sc
  enth.pm.y <- (integral2(fun = enth.post.y.sc, xmin = 0, xmax = 1, ymin = 0, ymax = 1)$Q)/enth.nc.sc
  
  # posterior means
  pm.mean.x <- inf.skpt.wt*skpt.pm.x + (1 - inf.skpt.wt)*enth.pm.x
  pm.mean.y <- inf.skpt.wt*skpt.pm.y + (1 - inf.skpt.wt)*enth.pm.y
  
  # coverage probability
  grid.x <- p.PC
  grid.y <- p.IP
  grid <- expand.grid(seq(0,1,by=0.01),seq(0,1,by=0.01))
  grid.index <- which.min(sqrt((grid$Var1-grid.x)^2+(grid$Var2-grid.y)^2))
  # evaluate "normalized density" at the grid
  grid.eval <- (inf.skpt.wt/skpt.nc.sc*skpt.post.sc(grid$Var1,grid$Var2) +
    (1-inf.skpt.wt)/enth.nc.sc*enth.post.sc(grid$Var1,grid$Var2))/
    sum(inf.skpt.wt/skpt.nc.sc*skpt.post.sc(grid$Var1,grid$Var2) +
        (1-inf.skpt.wt)/enth.nc.sc*enth.post.sc(grid$Var1,grid$Var2))
  #grid.eval[is.nan(grid.eval)]  <-  0 # 10-29-2019
  # find the cutoff for cred tail percentile (absolute value unnecessary)
  grid.index2 <- which.min(abs(cumsum(sort(grid.eval)) - cred.tail))
  # is grid point in credible interval?
  coverage <- grid.eval[grid.index] >= sort(grid.eval)[grid.index2]  
  
  return(cbind(pm.mean.x,pm.mean.y,coverage))
}

## Final Inference --------------------------------------------
monitoring <- function(index, fut.mix.prob, eff.mix.prob, inf.mix.prob){
  
  # efficacy and futility probabilities
  eff.prob.skpt <- integral2(skpt.post.sc, 
                             xmin = 0, 
                             xmax = 1 - delta.skpt,
                             ymin = function(x) x + delta.skpt,
                             ymax = 1)$Q
  eff.prob.enth <- integral2(enth.post.sc, 
                             xmin = 0, 
                             xmax = 1 - delta.skpt,
                             ymin = function(x) x + delta.skpt,
                             ymax = 1)$Q
  fut.prob.skpt <- 1 - integral2(skpt.post.sc, 
                                 xmin = 0,
                                 xmax = 1 - delta.intr,
                                 ymin = function(x) x + delta.intr,
                                 ymax = 1)$Q
  fut.prob.enth <- 1 - integral2(enth.post.sc, 
                                 xmin = 0, 
                                 xmax = 1-delta.intr,
                                 ymin = function(x) x + delta.intr,
                                 ymax = 1)$Q
  
  eff.prob <- (eff.skpt.wt/skpt.nc.sc)*eff.prob.skpt + ((1 - eff.skpt.wt)*enth.nc.sc)*eff.prob.enth
  inf.prob <- (inf.skpt.wt/skpt.nc.sc)*eff.prob.skpt + ((1 - inf.skpt.wt)*enth.nc.sc)*eff.prob.enth
  fut.prob <- (fut.skpt.wt/skpt.nc.sc)*fut.prob.skpt + ((1 - fut.skpt.wt)*enth.nc.sc)*fut.prob.enth
  
  return(cbind(eff.prob, inf.prob, fut.prob))
}