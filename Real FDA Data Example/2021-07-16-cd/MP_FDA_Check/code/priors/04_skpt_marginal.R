marginal <- function(){
  
  # grid search, check if final result is at end of interval
  alpha0.seq <- seq(1E-4, 2,length = 100)
  beta0.seq  <- seq(1,    6,length = 100)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      skpt.rd.alpha0 <- alpha0.seq[i]
      skpt.rd.beta0  <- beta0.seq[j]
      
      # closed form
      skpt.rd.prior <- function(x){
        dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
          (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
             pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))
      }
      
      # P(theta > theta_1) (called delta.enth for historical reasons)
      result1[i,j] <- integrate(skpt.rd.prior, 
                                lower = delta.enth, 
                                upper = 1)$value
      
      # P((theta_1+theta_0)/2 < theta < theta_1) (called delta.intr for historical reasons)
      result2[i,j] <- integrate(skpt.rd.prior, 
                                lower = delta.intr, 
                                upper = delta.enth)$value
    }
    
    # if (i%%10 == 0){print(paste0("Simulation ", i))}
    # 
    # # check normalizing constant
    # if (i%%10 == 0){print(paste0("Normalizing Constant ",
    #                              integrate(skpt.rd.prior,
    #                                        lower = -1,
    #                                        upper = 1)$value))}
  }
  
  result3 <- (result1 - q.outer)^2 + (result2 - q.inner)^2
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  
  print(paste0("Upper tail ", result1[i,j]))
  print(paste0("Upper half ", result2[i,j]))
  print(paste0("skpt.rd.alpha0 ", alpha0.seq[i]))
  print(paste0("skpt.rd.beta0  ", beta0.seq[j]))
  
  skpt.rd.alpha0 <- alpha0.seq[i]
  skpt.rd.beta0  <- beta0.seq[j]
  
  assign("skpt.rd.alpha0", skpt.rd.alpha0, envir = .GlobalEnv)
  assign("skpt.rd.beta0",  skpt.rd.beta0,  envir = .GlobalEnv)
}