##################################
### Evan Kwiatkowski, Feb 2020
###
### contains nested source("code_posteriors.R", local = TRUE)
###
##################################

# Posterior Mean & Coverage Probability  ---
pm_cp <- function(index){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)

  source("05_code_posteriors.R", local = TRUE)
  
  # posterior means
  skpt.pm.x <- (integrate_debug(skpt.post.x.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                integrate_debug(skpt.post.x.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
                skpt.nc.sc
  skpt.pm.y <- (integrate_debug(skpt.post.y.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                integrate_debug(skpt.post.y.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
                skpt.nc.sc
  enth.pm.x <- (integrate_debug(enth.post.x.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                integrate_debug(enth.post.x.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
                enth.nc.sc
  enth.pm.y <- (integrate_debug(enth.post.y.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                integrate_debug(enth.post.y.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
                enth.nc.sc
  
  # posterior means weighted
  pm.mean.x <- inf.skpt.wt*skpt.pm.x + (1 - inf.skpt.wt)*enth.pm.x
  pm.mean.y <- inf.skpt.wt*skpt.pm.y + (1 - inf.skpt.wt)*enth.pm.y
  
  ## coverage probability
  # create grid and find index of point representing true response probabilities
  grid.x <- p.IP - p.PC
  grid.y <- p.IP
  x.len <- 101
  y.len <- 101
  x <- seq(-1, 1, length = x.len)
  y <- seq(0,  1, length = y.len)
  # subset grid to trapizoid
  grid.outer <- expand.grid(x = x, y = y)
  grid <- grid.outer[(grid.outer$x >= 0 & (grid.outer$y < 1 - grid.outer$x)) | 
                     (grid.outer$x < 0  & (grid.outer$y > -grid.outer$x)), ]
  
  grid.index <- which.min(sqrt((grid$x - grid.x)^2+(grid$y - grid.y)^2))
  
  skpt.post.sc.grid <- NA
  enth.post.sc.grid <- NA
  for (i in 1:nrow(grid)){
    if (grid$x[i] >= 0 & (grid$y[i] < 1 - grid$x[i])) {
      skpt.post.sc.grid[i] <- skpt.post.sc.1(grid$x[i], grid$y[i])
      enth.post.sc.grid[i] <- enth.post.sc.1(grid$x[i], grid$y[i])
    } else if (grid$x[i] < 0 & (grid$y[i] > -grid$x[i])) {
      skpt.post.sc.grid[i] <- skpt.post.sc.2(grid$x[i], grid$y[i])
      enth.post.sc.grid[i] <- enth.post.sc.2(grid$x[i], grid$y[i])
    } else {
      skpt.post.sc.grid[i] <- 0
      enth.post.sc.grid[i] <- 0
    }
  }
  
  # Note: grid.eval crashes if any cell count is 0
  skpt.post.sc.grid[is.nan(skpt.post.sc.grid)]  <-  0 # 10-29-2019, only if necessary
  enth.post.sc.grid[is.nan(enth.post.sc.grid)]  <-  0 # 10-29-2019, only if necessary
  
  # evaluate "normalized density" at the grid
  grid.eval <- (inf.skpt.wt*skpt.post.sc.grid/skpt.nc.sc + (1 - inf.skpt.wt)*enth.post.sc.grid/enth.nc.sc)/
            sum(inf.skpt.wt*skpt.post.sc.grid/skpt.nc.sc + (1 - inf.skpt.wt)*enth.post.sc.grid/enth.nc.sc)

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
  fut.prob.skpt <- integrate_debug(skpt.post.sc.1, xmin = delta.intr, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
  
  # futility probability using enthusiastic prior (scaled)
  # integrate joint prior from [delta.intr , 1], stop if less than epsilon
  # alternatively, stop if 1 - (fut.prob.skpt/skpt.nc.sc) is greater than 1 - epsilon
  fut.prob.enth <- integrate_debug(enth.post.sc.1, xmin = delta.intr, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
  
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
  
  return(data.frame(cbind(eff.prob, inf.prob, fut.prob, eff.mix.prob, risk.diff.mle, skpt.psi, enth.psi, ni.psi)))
}