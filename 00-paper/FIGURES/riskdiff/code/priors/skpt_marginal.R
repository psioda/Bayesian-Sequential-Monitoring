marginal <- function(){
  
  alpha0.seq <- seq(1E-4, 2,length = 100)
  beta0.seq  <- seq(1,    6,length = 100)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      skpt.rd.alpha0 <- alpha0.seq[i]
      skpt.rd.beta0  <- beta0.seq[j]
      
      skpt.rd.prior <- function(x){
        exp(-(abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0)*
          skpt.rd.beta0/(2*skpt.rd.alpha0*gamma(1/skpt.rd.beta0))/
          (pgnorm(q = 1,  mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
          pgnorm(q = -1, mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))
      }
      
      result1[i,j] <- integrate(skpt.rd.prior, lower = delta.enth, upper = 1)$value
      result2[i,j] <- integrate(skpt.rd.prior, lower = delta.intr, upper = delta.enth)$value
    }
    if (i%%10 == 0){print(paste0("Simulation ", i))}
    if (i%%10 == 0){print(integrate(skpt.rd.prior, lower = -1, upper = 1)$value)}
  }
  
  result3 <- (result1 - q.outer)^2 + (result2 - q.inner)^2
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  
  print(result1[i,j])
  print(result2[i,j])
  
  skpt.rd.alpha0 <- alpha0.seq[i]
  skpt.rd.beta0 <- beta0.seq[j]
  
  assign("skpt.rd.alpha0",skpt.rd.alpha0,envir = .GlobalEnv)
  assign("skpt.rd.beta0",skpt.rd.beta0,envir = .GlobalEnv)
}