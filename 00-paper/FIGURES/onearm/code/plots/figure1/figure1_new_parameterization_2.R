rm(list = ls())
library(gnorm)
## equivalence between dnorm and dgnorm
mu0.enth     <- runif(1)
sigma0.enth  <- runif(1)
lambda0.enth <- 2
dnorm(mu0.enth,  mean = mu0.enth, sd = sqrt(sigma0.enth))
dgnorm(mu0.enth, mu = mu0.enth, alpha = sqrt(2*sigma0.enth), beta = lambda0.enth)

rm(list = ls())
## standard normal
a <- dnorm(0)
b <- pnorm(qnorm(0.975))

enth_prior_default <- function(){
  mu0.enth     <- p.enth
  sigma0.seq   <- seq(.01, 0.5, by=0.0001)
  lambda0.enth <- 2
  result       <- NA
  
  for (i in 1:length(sigma0.seq)){
    sigma0.enth <- sigma0.seq[i]
    result[i] <- pgnorm(p.skpt, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)
  }
  
  i <- which(abs(result-tail.enth)==min(abs(result-tail.enth)))
  sigma0.enth <- sigma0.seq[i]
  prior.enth <- function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
  }
  nc.enth <- integrate(prior.enth, lower=-Inf, upper=Inf)[[1]]
  prior.nc.enth <- function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
  }
  print(paste0("mu: ", mu0.enth, ",  sigma: ", sigma0.enth, ",  lambda: ", lambda0.enth))
  print(paste0("Tail area: ", result[i]))
  print(paste0("Half-width area: ", 
               integrate(prior.nc.enth, lower=p.skpt, upper=p.intr)[[1]]))
  assign("mu0.enth", mu0.enth, envir = .GlobalEnv)
  assign("sigma0.enth", sigma0.enth, envir = .GlobalEnv)
  assign("lambda0.enth", lambda0.enth, envir = .GlobalEnv)
  return(prior.nc.enth)
}

p.skpt    <- 0.40     # response rate for skeptic,  enthusiast,  futility
p.enth    <- 0.67
p.intr    <- (p.skpt+p.enth)/2
tail.skpt <- 0.025  # tail probabilities for priors
tail.enth <- 0.025
epsilon   <- 0
prior.nc.enth <- enth_prior_default()

sigma0.enth.normal <- sigma0.enth

enth_prior_custom <- function(scale, method){
  
  # generate 2 probability conditions over the grid
  mu0.enth    <- p.enth
  sigma0.seq  <- seq(.01, 1, by=0.001)
  lambda0.seq <- seq(1, 7, by=0.01)
  result1 <- matrix(NA, nrow=length(sigma0.seq), ncol=length(lambda0.seq))
  result2 <- matrix(NA, nrow=length(sigma0.seq), ncol=length(lambda0.seq))
  result2b <- matrix(NA, nrow=length(sigma0.seq), ncol=length(lambda0.seq))
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      sigma0.enth   <- sigma0.seq[i]
      lambda0.enth  <- lambda0.seq[j]
      result1[i, j] <- pgnorm(p.skpt, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)
      if (method == 1){result2[i, j] <- pgnorm(p.intr, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth) - 
        pgnorm(p.skpt, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)}
      if (method == 2){result2[i, j] <- dgnorm(0, mu = 0, alpha = sigma0.enth, beta = lambda0.enth)}
    }
  }
  
  # find closest match and assign variables
  if (method == 1){result3 <- abs(result1 - tail.enth) + abs(result2 - (pnorm(qnorm(tail.enth)/2) - tail.enth)*scale)}
  if (method == 2){result3 <- abs(result1 - tail.enth) + abs(result2 - dnorm(0, mean = 0, sd = sigma0.enth.normal/sqrt(2))*scale)}
  index   <- which(result3 == min(result3),  arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  sigma0.enth <- sigma0.seq[i]
  lambda0.enth <- lambda0.seq[j]
  prior.nc.enth <- function(x){
    dgnorm(x, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)
  }
  print(paste0("mu: ", mu0.enth, ",  sigma: ", sigma0.enth, ",  lambda: ", lambda0.enth))
  print(paste0("Tail area: ", result1[i, j]))
  print(paste0("Tail area 2: ", pgnorm(p.intr, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth) - 
                 pgnorm(p.skpt, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)))
  if (method == 1){print(paste0("Tail area 2 condition: ", (pnorm(qnorm(tail.enth)/2) - tail.enth)*scale))}
  if (method == 2){print(paste0("Density condition: ", dnorm(0, mean = 0, sd = sigma0.enth.normal/sqrt(2))*scale))}
  print(paste0("Density at mode: ", dgnorm(0, mu = 0, alpha = sigma0.enth, beta = lambda0.enth)))
  assign("mu0.enth", mu0.enth, envir = .GlobalEnv)
  assign("sigma0.enth", sigma0.enth, envir = .GlobalEnv)
  assign("lambda0.enth", lambda0.enth, envir = .GlobalEnv)
  assign("tail.enth.actual", result1[i, j], envir = .GlobalEnv)
  return(prior.nc.enth)
}
#scale  <-  0.6844146 # for method 2 to match scale(1.5) of method 1
#scale  <- 1.433264 # for method 2 to match scale(0.75) of method 1
scale  <- 1.5
method <- 1
prior.nc.enth <- enth_prior_custom(scale=scale, method=method)

# check density at intermediate value
dgnorm(p.intr, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)

# compare to old system
pgnorm(p.intr, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth) - 
  pgnorm(p.skpt, mu = mu0.enth, alpha = sigma0.enth, beta = lambda0.enth)

# test density values
dgnorm(0, mu = 0, alpha = sigma0.enth, beta = lambda0.enth)
dnorm(0, mean = 0, sd = sigma0.enth/sqrt(2))


xmin  <-  p.enth - 0.5
xmax  <-  p.enth + 0.5
ymax  <-  7
x <- seq(xmin, 
         xmax, 
         by=0.005)
plot(x, prior.nc.enth(x), type="l", 
     xlab="", 
     ylab="", 
     main="", 
     xaxt="n", 
     yaxt="n", 
     xlim=c(xmin, xmax), 
     ylim=c(0, ymax)) # 20-01-02
axis(2, at = seq(0, ymax, by = 1))
axis(1, at=c(p.enth, p.skpt, (p.enth+p.skpt)/2), 
     labels=c(as.expression(bquote(theta[1])), 
              as.expression(bquote(theta[0])), 
              as.expression(bquote((theta[0]+theta[1])/2))))
title(ylab="Density Value",  line=2)
title(xlab="Response Probability", line=2)
#title(xlab="Density Value", line=2)
polygon(c(x[x<=p.skpt], p.skpt), 
        c(prior.nc.enth(x)[x<=p.skpt], 0), col="black")
polygon(c(p.skpt, x[x>=p.skpt & x<=p.intr], p.intr), 
        c(0, prior.nc.enth(x)[x>=p.skpt & x<=p.intr], 0), col="lightgrey")
segments(x0=p.enth, y0=0, y1=prior.nc.enth(p.enth))
legend("top", 
       legend= c(
         as.expression(bquote(mode(theta) == theta[1])), 
         as.expression(bquote(P(theta< theta[0])==.(tail.enth))), 
         as.expression(bquote(P(theta %in% (theta[0]*", "*(theta[0]+theta[1])/2)==.(round((pnorm(qnorm(tail.enth)/2)-tail.enth)*scale, 3)))))#, 
         #as.expression(bquote(GN(mu==theta[1], alpha==.(sigma0.enth), beta==.(lambda0.enth))))
       ))

