fcn_prior_placebo <- function(){
  
  sigma0.seq <- seq(1E-2, 1, length = 10000)
  lambda0.seq <- 5
  
  result1 <- matrix(NA, nrow = length(sigma0.seq), ncol = length(lambda0.seq))
  result2 <- matrix(NA, nrow = length(sigma0.seq), ncol = length(lambda0.seq))
  
  q.outer <- 0.5  # window around mu +/- delta.enth
  q.inner <- 0.25 # window around mu +/- delta.intr
  
  for (i in 1:length(sigma0.seq)){
    if (i%%1000 == 0){print(paste0("Simulation ", i))}
    
    for (j in 1:length(lambda0.seq)){
      
      placebo.sigma0 <- sigma0.seq[i]
      placebo.lambda0 <- lambda0.seq[j]
      
      placebo.prior <- function(x){
        exp(-(abs(x-mu)/placebo.sigma0)^placebo.lambda0)
      }
      
      placebo.nc <- integrate(placebo.prior, 
                      lower=0, 
                      upper=1)$value

      result1[i,j] <- (integrate(placebo.prior,
                                lower = mu - delta.enth,
                                upper = mu + delta.enth)$value)/placebo.nc
      result2[i,j] <- (integrate(placebo.prior,
                                lower = mu - delta.intr,
                                upper = mu + delta.intr)$value)/placebo.nc
    }
  }
  
  result3 <- abs(result1 - q.outer) + abs(result2 - q.inner)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  
  print(result1[i,j])
  print(result2[i,j])
  
  placebo.sigma0 <- sigma0.seq[i]
  placebo.lambda0 <- lambda0.seq[j]
  
  assign("placebo.sigma0",placebo.sigma0,envir = .GlobalEnv)
  assign("placebo.lambda0",placebo.lambda0,envir = .GlobalEnv)
  
}