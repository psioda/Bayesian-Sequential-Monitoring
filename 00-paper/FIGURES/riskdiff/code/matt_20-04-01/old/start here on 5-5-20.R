rm(list=ls())

# Evaluates prior-data conflict metric 
# Evans2006, pg897, ex2
evans <- function(n, t0, a, b){
  
  post.nc      <- NA
  for (x in 1:(n+1)){ # has to start binomial looping at 0
    post.kernel   <- function(p) dbinom(x-1, n, p)*dbeta(p, a, b)
    post.nc[x]    <- integrate(post.kernel, 0, 1)[[1]]
  }
  
  sum(post.nc[post.nc <= post.nc[t0 + 1]])
  
}


evans(10,5,4,8)

rm(list=ls())
library(pracma)

evans2d <- function(n1, x1, a1, b1, 
                    n2, x2, a2, b2){

    post.nc       <- matrix(NA, nrow = n1 + 1, ncol = n2 + 1)
  

  for (i in 1:(n1+1)){
    for (j in 1:(n2+1)){
      #post.nc <- NA
      post.kernel   <- function(p1, p2) dbinom(i-1, n1, p1)*
                                        dbinom(j-1, n2, p2)*
                                        dbeta(p1, a1, b1)*
                                        dbeta(p2, a2, b2)
      #post.nc <- integral2(post.kernel, 1E-1, 1, 1E-1, 1, singular = TRUE)[[1]]
      post.nc[i, j] <- integral2(post.kernel, 0, 1, 0, 1, singular = TRUE)[[1]]
    }
  }
  sum(post.nc[post.nc <= post.nc[x1 + 1, x2 + 1]])
}
evans2d_log <- function(n1, x1, a1, b1, 
                        n2, x2, a2, b2){
  
  post.nc       <- matrix(NA, nrow = n1 + 1, ncol = n2 + 1)
  
  
  for (i in 1:(n1+1)){
    for (j in 1:(n2+1)){
      #post.nc <- NA
      post.kernel   <- function(p1, p2) 
        dbinom(i-1, n1, p1, log = TRUE) +
        dbinom(j-1, n2, p2, log = TRUE) +
        dbeta(p1, a1, b1, log = TRUE) +
        dbeta(p2, a2, b2, log = TRUE)
      #post.nc <- integral2(post.kernel, 1E-1, 1, 1E-1, 1, singular = TRUE)[[1]]
      post.nc[i, j] <- integral2(post.kernel, 0, 1, 0, 1, singular = TRUE)[[1]]
    }
  }
  sum(post.nc[post.nc <= post.nc[x1 + 1, x2 + 1]])
}


# seems to be working
evans2d(10,3,4,8,20,12,8,4)
evans2d(10,8,4,8,20,4,8,4)
# next use post.nc to verify a value...




# goal: work up to showing that mixing weight (things from nc.blah blah) are the same
# will need to see if R understands using Rbinom with (x+y) as response probability, then integrating

# is not the end of the world if it doesn't work

# or could skip this step for now ... 


# big picture anything that shows the 2d prior predictive distribution and prior data conflict calculation are working

evans <- function(n, t0, a, b){
  
  n<-10
  a<-4
  b<-8
  t0<-5
  
  prior.kernel <- function(p) p^(a-1)*(1-p)^(b-1)
  prior.nc     <- integrate(prior.kernel, 0, 1)[[1]]
  
  
  post.nc      <- NA
  for (x in 1:(n+1)){ # has to start binomial looping at 0
    post.kernel   <- function(p) choose(n,(x-1))*p^((x-1)+a-1)*(1-p)^(n-(x-1)+b-1)/prior.nc
    post.nc[x]    <- integrate(post.kernel, 0, 1)[[1]]
  }
  
  sum(post.nc[post.nc <= post.nc[t0 + 1]])
  
}


# scaled (un-normalized) posterior density
skpt.post.sc.1 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
            pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}


test <- function(p1,p2){
  dbinom(x, n, p1+p2)
}
n<-10
x<-5
library("pracma")
integrate_debug(test,0,1,0,0.5)

