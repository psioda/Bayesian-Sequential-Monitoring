##################################
### Evan Kwiatkowski, Feb 2020
##################################

# Posterior Mean & Coverage Probability  ---
pm_cp <- function(index){
  
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
  grid <- expand.grid(seq(0, 1,by = 0.01), seq(0,1,by = 0.01))
  grid.index <- which.min(sqrt((grid$Var1 - grid.x)^2+(grid$Var2 - grid.y)^2))
  
  skpt.post.sc.grid <- skpt.post.sc(grid$Var1, grid$Var2)
  enth.post.sc.grid <- enth.post.sc(grid$Var1, grid$Var2)
  # Note: grid.eval crashes if any cell count is 0
  skpt.post.sc.grid[is.nan(skpt.post.sc.grid)]  <-  0 # 10-29-2019, only if necessary
  enth.post.sc.grid[is.nan(enth.post.sc.grid)]  <-  0 # 10-29-2019, only if necessary
  
  # evaluate "normalized density" at the grid
  grid.eval <- (inf.skpt.wt*skpt.post.sc.grid/skpt.nc.sc +
    (1 - inf.skpt.wt)*enth.post.sc.grid/enth.nc.sc)/
    sum(inf.skpt.wt*skpt.post.sc.grid/skpt.nc.sc +
        (1 - inf.skpt.wt)*enth.post.sc.grid/enth.nc.sc)

  # find the cutoff for cred tail percentile (absolute value unnecessary)
  grid.index2 <- which.min(abs(cumsum(sort(grid.eval)) - cred.tail))
  
  # is grid point in credible interval?
  coverage <- (grid.eval[grid.index] >= sort(grid.eval)[grid.index2]) 
  
  mle.IP <- y1.IP/sum(y0.IP, y1.IP)
  mle.PC <- y1.PC/sum(y0.PC, y1.PC)
  
  return(data.frame(cbind(pm.mean.x, pm.mean.y, coverage, mle.PC, mle.IP)))
}

# Monitoring ---
monitoring <- function(index){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)
  
  source("code_posteriors.R", local = TRUE)
  
  # efficacy probability using skeptical prior
  eff.prob.skpt <- integrate_debug(skpt.post.sc, 
                             xmin = 0, 
                             xmax = 1 - delta.skpt,
                             ymin = function(x) x + delta.skpt,
                             ymax = 1)
  # efficacy probability using enthusiastic prior
  eff.prob.enth <- integrate_debug(enth.post.sc, 
                             xmin = 0, 
                             xmax = 1 - delta.skpt,
                             ymin = function(x) x + delta.skpt,
                             ymax = 1)
  # futility probability using skeptical prior
  fut.prob.skpt <- integrate_debug(skpt.post.sc, 
                                 xmin = 0,
                                 xmax = 1 - delta.intr,
                                 ymin = function(x) x + delta.intr,
                                 ymax = 1)
  # futility probability using enthusiastic prior
  fut.prob.enth <- integrate_debug(enth.post.sc, 
                                 xmin = 0, 
                                 xmax = 1 - delta.intr,
                                 ymin = function(x) x + delta.intr,
                                 ymax = 1)
  
  eff.prob <- eff.skpt.wt*(eff.prob.skpt/skpt.nc.sc) + (1 - eff.skpt.wt)*(eff.prob.enth/enth.nc.sc)
  inf.prob <- inf.skpt.wt*(eff.prob.skpt/skpt.nc.sc) + (1 - inf.skpt.wt)*(eff.prob.enth/enth.nc.sc)
  fut.prob <- fut.skpt.wt*((1 - fut.prob.enth/enth.nc.sc) + (1 - fut.skpt.wt)*(1 - fut.prob.enth/enth.nc.sc))
  # bug caught 2/1/20, (1 - fut.prob.skpt)) must be computed before dividing by normalizing constant
  # bug caught 2/2/20, 1 - fut.prob.enth/enth.nc.sc NOT (1 - fut.prob.enth)/enth.nc.sc
  
  return(data.frame(cbind(eff.prob, inf.prob, fut.prob)))
}
