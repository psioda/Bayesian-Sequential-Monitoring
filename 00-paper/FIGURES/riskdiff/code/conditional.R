skpt_tail_area <- function(){
  
  alpha0.seq <- seq(1E-1,  1, length = 10)
  beta0.seq  <- seq(1E-1,    1, length = 10)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      skpt.alpha0 <- alpha0.seq[i]
      skpt.beta0  <- beta0.seq[j]
      
      # normalized prior density
      skpt.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
        exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)/(2*placebo.alpha0*gamma(1/placebo.beta0)/placebo.beta0)*
          exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
          (pgnorm(q = 1,  mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0) -
             pgnorm(q = -1, mu = delta.skpt, alpha = placebo.alpha0, beta = placebo.beta0))/
          (pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
             pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
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
      c1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu - delta.enth, ymax = mu + delta.enth)
      c2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu - delta.enth, ymax = mu + delta.enth)
      result1[i,j] <- c1 + c2
      d1 <- integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = mu - delta.intr, ymax = mu + delta.intr)
      d2 <- integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = mu - delta.intr, ymax = mu + delta.intr)
      result2[i,j] <- d1 + d2
    }
    if (i%%10 == 0){print(paste0("Simulation ", i))}
    if (i%%10 == 0){print(
      integrate_debug(skpt.prior.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
      integrate_debug(skpt.prior.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
    )}
    }
  
  result3 <- sqrt((result1 - q.outer)^2)  + sqrt((result2 - q.inner)^2)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  skpt.alpha0 <- alpha0.seq[i]
  skpt.beta0 <- beta0.seq[j]
  print(result1[i,j])
  print(result2[i,j])
  
  assign("skpt.alpha0", skpt.alpha0, envir = .GlobalEnv)
  assign("skpt.beta0", skpt.beta0, envir = .GlobalEnv)
  
}