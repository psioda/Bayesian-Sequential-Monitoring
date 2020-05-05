rm(list=ls())

n <- 10
x <- 5
a <- 4
b <- 8
m <- a/(a+b)
s <- (a+b)

prior.kernel  <- function(p) p^(a-1)*(1-p)^(b-1)
prior.nc      <- integrate(prior.kernel, 0, 1)[[1]]

post.kernel   <- function(p) choose(n,x)*p^(x+a-1)*(1-p)^(n-x+b-1)/prior.nc
post.nc       <- integrate(post.kernel, 0, 1)[[1]]
post.density  <- function(p) post.kernel(p)/post.nc

# Verify that posterior distribution can be calculated computationally
p.seq        <- seq(0,1,length=5)
norm(post.density(p.seq)-dbeta(p.seq,x+a,n-x+b),type="2")

# Verify that prior predictive distribution can be calculated computationally
library(rmutil)
norm(post.nc-dbetabinom(x,n,m,s),type="2")


rm(list=ls())

# Evaluates prior-data conflict metric 
# Uses closed-form prior predictive values
# Evans2006, pg897, ex2
evans <- function(n, t0, a, b){
  
  m <- a/(a+b)
  s <- (a+b)
  
  
  
  sum(dbetabinom(0:n, size = n, m = m, s = s)[
    round(dbetabinom(0:n, size = n, m = m, s = s),5) <= 
    round(dbetabinom(t0, size = n, m = m, s = s),5)])
}

n <- 13
a.skpt <- 4
b.skpt <- 8
a.enth <- 8
b.enth <- 4

psi.skpt <- NA
psi.enth <- NA

for (i in 1:(n+1)){
  psi.skpt[i] <- evans(n, i-1, a.skpt, b.skpt)
  psi.enth[i] <- evans(n, i-1, a.enth, b.enth)
}

psi.enth-psi.skpt

plot(0:n,psi.skpt,type='l')
lines(0:n,psi.enth,type='l')

rm(list=ls())

y <- NA
n <- 100
a <- 5
b <- 20
for (i in 1:n+1){
  y[i] <- evans(n = n, i-1,a=a,b=b)
}
