ni_tail_area <- function(){
  
  # grid search, check if final result is at end of interval
  alpha0.seq <- seq(0.125, 0.5,  length = 20)
  beta0.seq  <- seq(1.8,   2.2,  length = 20)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      ni.alpha0 <- alpha0.seq[i]
      ni.beta0  <- beta0.seq[j]
      
      # normalized prior density
      ni.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
        dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
          (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
             pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
          dgnorm(y,          mu,        ni.alpha0,    ni.beta0)/
          (pgnorm(q = 1 - x, mu,        ni.alpha0,    ni.beta0) -
             pgnorm(q = 0,   mu,        ni.alpha0,    ni.beta0))
      }
      ni.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
        dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
          (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
             pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
          dgnorm(y,         mu,         ni.alpha0,    ni.beta0)/
          (pgnorm(q = 1,    mu,         ni.alpha0,    ni.beta0) -
             pgnorm(q = -x, mu,         ni.alpha0,    ni.beta0))
      }

      # tail probabilities
      c1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
      c2 <- integrate_debug(ni.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
      c3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
      result1[i,j] <- c1 + c2 + c3
      d1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
      d2 <- integrate_debug(ni.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
      d3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
      result2[i,j] <- d1 + d2 + d3 - c1 - c2 - c3
    }
    
    # if (i%%10 == 0){print(paste0("Simulation ", i))}
    # 
    # # check normalizing constant
    # if (i%%10 == 0){print(paste0("Normalizing Constant ", 
    #   integrate_debug(ni.prior.1, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x) +
    #   integrate_debug(ni.prior.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
    # ))}
    }
  
  result3 <- sqrt((result1 - q.outer)^2)  + sqrt((result2 - q.inner)^2)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  
  ni.alpha0 <- alpha0.seq[i]
  ni.beta0  <- beta0.seq[j]
  
  print(paste0("Upper tail ", result1[i,j]))
  print(paste0("Upper half ", result2[i,j]))
  print(paste0("ni.alpha0 ", alpha0.seq[i]))
  print(paste0("ni.beta0  ", beta0.seq[j]))
  
  assign("ni.alpha0", ni.alpha0, envir = .GlobalEnv)
  assign("ni.beta0", ni.beta0,   envir = .GlobalEnv)
}