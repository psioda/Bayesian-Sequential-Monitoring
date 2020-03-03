rm(list = ls())

library(pracma)

scale <- c(0.75, 1, 1.5)
scale.a <- NA
scale.b <- NA

for (k in 1:length(scale)){
delta.enth <- 1.959964
delta.skpt <- 0
delta.intr <- (delta.skpt+delta.enth)/2

q.outer    <- 0.025               # y > x + delta.enth
q.inner    <- (1 - 0.8364525 - 0.025)*scale[k]   # y > x + delta.intr (major typo caught!)

marginal <- function(){
  
  alpha0.seq <- seq(0.2, 6, length = 100)
  beta0.seq  <- seq(0.5, 6, length = 100)
  
  result1 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  result2 <- matrix(NA, nrow = length(alpha0.seq), ncol = length(beta0.seq))
  
  for (i in 1:length(alpha0.seq)){
    for (j in 1:length(beta0.seq)){
      
      placebo.alpha0 <- alpha0.seq[i]
      placebo.beta0  <- beta0.seq[j]
      
      placebo.prior <- function(x){
        exp(-(abs(x - delta.skpt)/placebo.alpha0)^placebo.beta0)*
          placebo.beta0/(2*placebo.alpha0*gamma(1/placebo.beta0))
      }
      
      result1[i,j] <- integrate(placebo.prior, lower = delta.enth, upper = Inf)$value
      result2[i,j] <- integrate(placebo.prior, lower = delta.intr, upper = delta.enth)$value
    }
    if (i%%100 == 0){print(paste0("Simulation ", i))}
    if (i%%100 == 0){print(integrate(placebo.prior, lower = -Inf, upper = Inf)$value)}
  }
  
  result3 <- (result1 - q.outer)^2 + (result2 - q.inner)^2
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  
  print(result1[i,j])
  print(result2[i,j])
  
  placebo.alpha0 <- alpha0.seq[i]
  placebo.beta0  <- beta0.seq[j]
  
  assign("placebo.alpha0",placebo.alpha0,envir = .GlobalEnv)
  assign("placebo.beta0", placebo.beta0, envir = .GlobalEnv)
  
}
marginal()

f1 <- function(a){   # q.outer
  placebo.prior <- function(x){
    exp(-(abs(x - delta.skpt)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))
  }
  integrate(placebo.prior, lower = delta.enth, upper = Inf)$value
}
f2 <- function(a){   # q.inner
  placebo.prior <- function(x){
    exp(-(abs(x - delta.skpt)/a[1])^a[2])*
      a[2]/(2*a[1]*gamma(1/a[2]))
  }
  integrate(placebo.prior, lower = delta.intr, upper = delta.enth)$value
}

start <- c(placebo.alpha0,placebo.beta0)
#start <- c (0.1, 2)
f1(start)
f2(start)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

nlm.fit <- nlminb(start     = start, 
                  objective = fn, 
                  lower     = c(0,0), 
                  upper     = c(Inf, Inf))
a2 <- nlm.fit$par
a2
f1(a2)
f2(a2)
scale.a[k] <- a2[1]
scale.b[k] <- a2[2]

}

x<-seq(-3,3,by=0.01)
k <- 1
placebo.prior <- function(x){
  exp(-(abs(x - delta.skpt)/scale.a[k])^scale.b[k])*
    scale.b[k]/(2*scale.a[k]*gamma(1/scale.b[k]))
}
plot(x,placebo.prior(x),type="l",
     xlab="",
     ylab="",
     main="",
     xaxt="n",
     yaxt="n")

for (k in 1:length(scale)){
placebo.prior <- function(x){
  exp(-(abs(x - delta.skpt)/scale.a[k])^scale.b[k])*
    scale.b[k]/(2*scale.a[k]*gamma(1/scale.b[k]))
}
lines(x,placebo.prior(x))
}

axis(1,at=c(delta.enth,delta.skpt),
     labels=c(as.expression(bquote(theta[1])),as.expression(bquote(theta[0]))))
title(ylab="Density Value", line=1, #cex.lab=1.2
)
title(xlab="Response Probability",line=2)

legend('topleft',
       legend= c(as.expression(bquote(mu == theta[0])),
                 #as.expression(bquote(alpha == .(sigma0.enth))),
                 #as.expression(bquote(beta == .(lambda0.enth))),
                 as.expression(bquote(P(theta> theta[1])==.(q.outer)))))
