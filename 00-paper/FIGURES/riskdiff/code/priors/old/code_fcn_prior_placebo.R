fcn_prior_placebo <- function(){
  
  alpha0.seq <- seq(1E-4, 10, length = 1000)
  beta0.seq  <- seq(1E-1, 4,  length = 10)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      placebo.alpha0 <- alpha0.seq[i]
      placebo.beta0  <- beta0.seq[j]
      
      placebo.prior <- function(x){
        exp(-(abs(x-mu)/placebo.alpha0)^placebo.beta0)*
          placebo.beta0/(2*placebo.alpha0*gamma(1/placebo.beta0))/
          (pgnorm(q = 1, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0) -
           pgnorm(q = 0, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0))
      }

      result1[i,j] <- (integrate(placebo.prior,
                                lower = mu - delta.enth,
                                upper = mu + delta.enth)$value)
      result2[i,j] <- (integrate(placebo.prior,
                                lower = mu - delta.intr,
                                upper = mu + delta.intr)$value)
    }
    if (i%%100 == 0){print(paste0("Simulation ", i))}
    if (i%%100 == 0){print(integrate(placebo.prior, lower = 0, upper = 1)$value)}
  }
  
  result3 <- (result1 - q.outer)^2 + (result2 - q.inner)^2
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  
  print(result1[i,j])
  print(result2[i,j])
  
  placebo.alpha0 <- alpha0.seq[i]
  placebo.beta0 <- beta0.seq[j]
  
  assign("placebo.alpha0",placebo.alpha0,envir = .GlobalEnv)
  assign("placebo.beta0",placebo.beta0,envir = .GlobalEnv)
  
}