enth_tail_area <- function(){
  
  alpha0.seq <- seq(1E-3, 1, length = 1000)
  enth.beta0 <- 2
  result1 <- matrix(NA, ncol=length(alpha0.seq))

  q.outer <- 0.025  # y > x + delta.enth

  for (i in 1:length(alpha0.seq)){
    if (i%%100 == 0){print(paste0("Simulation ", i))}
    
      enth.alpha0 <- alpha0.seq[i]

      # (un-normalized) prior density
      enth.prior <- function(x, y){
        exp(
          - (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
            (abs((y - x) - delta.enth)/enth.alpha0)^enth.beta0
        )
      }
      
      # normalizing constant for prior density
      enth.prior.nc <- integrate_debug(fun = enth.prior,
                                       xmin = 0,
                                       xmax = 1,
                                       ymin = 0,
                                       ymax = 1)
      
      result1[i] <- integrate_debug(fun = enth.prior,
                                      xmin = 0,
                                      xmax = 1 - delta.skpt,
                                      ymin = function(x) x + delta.skpt,
                                      ymax = 1)/enth.prior.nc
  }

  result3 <- abs(result1 - (1 - q.outer))
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[2] # Jan 30 correction
  enth.alpha0 <- alpha0.seq[i]

  print(1 - result1[i])

  assign("enth.alpha0", enth.alpha0, envir = .GlobalEnv)
  assign("enth.beta0", enth.beta0, envir = .GlobalEnv)
}