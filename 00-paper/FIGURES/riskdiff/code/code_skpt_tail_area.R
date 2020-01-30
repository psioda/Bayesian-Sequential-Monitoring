skpt_tail_area <- function(){
  
  sigma0.seq <- seq(1E-2, 1, length = 1000)
  lambda0.seq <- 1.5
  result1 <- matrix(NA, nrow = length(sigma0.seq), ncol = length(lambda0.seq))
  result2 <- matrix(NA, nrow = length(sigma0.seq), ncol = length(lambda0.seq))
  
  q.outer <- 0.025  # y > x + delta.enth
  #q.inner <- 0.10  # y > x + delta.intr ( major typo caught! )
  
  for (i in 1:length(sigma0.seq)){
    if (i%%100 == 0){print(paste0("Simulation ", i))}
    
    for (j in 1:length(lambda0.seq)){
      
      skpt.sigma0 <- sigma0.seq[i]
      skpt.lambda0 <- lambda0.seq[j]
      
      # (un-normalized) prior density
      skpt.prior <- function(x, y){
        exp(
          - (abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
            (abs((y - x) - delta.skpt)/skpt.sigma0)^skpt.lambda0
        )
      }
      # normalizing constant for prior density
      skpt.prior.nc <- integrate_debug(fun = skpt.prior,
                                    xmin = 0,
                                    xmax = 1,
                                    ymin = 0,
                                    ymax = 1)

      # tail probabilities
      result1[i,j] <- integrate_debug(fun = skpt.prior,
                                     xmin = 0,
                                     xmax = 1 - delta.enth,
                                     ymin = function(x) x + delta.enth,
                                     ymax = 1)/skpt.prior.nc
      # result2[i,j] <- integrate_debug(fun = skpt.prior,
      #                                xmin = 0,
      #                                xmax= 1 - delta.intr,
      #                                ymin = function(x) x + delta.intr,
      #                                ymax = 1)/skpt.prior.nc
    }
  }
  
  result3 <- sqrt((result1 - q.outer)^2) # + sqrt((result2 - q.inner)^2)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  skpt.sigma0 <- sigma0.seq[i]
  skpt.lambda0 <- lambda0.seq[j]
  print(result1[i,j])
  # print(result2[i,j])
  
  assign("skpt.sigma0", skpt.sigma0, envir = .GlobalEnv)
  assign("skpt.lambda0", skpt.lambda0, envir = .GlobalEnv)
  
}