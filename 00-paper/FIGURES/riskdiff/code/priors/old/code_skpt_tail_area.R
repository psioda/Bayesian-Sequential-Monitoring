skpt_tail_area <- function(){
  
  alpha0.seq <- seq(1E-6, 2, length = 100)
  beta0.seq  <- seq(0.1, 10, length = 10)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      skpt.alpha0 <- alpha0.seq[i]
      skpt.beta0  <- beta0.seq[j]
      
      # normalized prior density
      skpt.prior <- function(x, y){
        exp(
          - (abs(x - mu)/placebo.alpha0)^placebo.beta0 - 
            (abs(y - (x + delta.skpt))/skpt.alpha0)^skpt.beta0)*
          placebo.beta0/(2*placebo.alpha0*gamma(1/placebo.beta0))*
          skpt.beta0/(2*skpt.alpha0*gamma(1/skpt.beta0))/
          (pgnorm(q = 1, mu = x + delta.skpt, alpha = skpt.alpha0, beta = skpt.beta0) -
           pgnorm(q = 0, mu = x + delta.skpt, alpha = skpt.alpha0, beta = skpt.beta0))/
          (pgnorm(q = 1, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0) -
           pgnorm(q = 0, mu = mu, alpha = placebo.alpha0, beta = placebo.beta0))
      }

      # tail probabilities
      result1[i,j] <- integrate_debug(fun = skpt.prior,
                                     xmin = 0,
                                     xmax = 1 - delta.enth,
                                     ymin = function(x) x + delta.enth,
                                     ymax = 1)
      result2[i,j] <- integrate_debug(fun = skpt.prior,
                                      xmin = 0,
                                      xmax= 1 - delta.intr,
                                      ymin = function(x) x + delta.intr,
                                      ymax = 1)
    }
    if (i%%10 == 0){print(paste0("Simulation ", i))}
    if (i%%10 == 0){print(integrate_debug(skpt.prior, 0, 1, 0 ,1))}
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