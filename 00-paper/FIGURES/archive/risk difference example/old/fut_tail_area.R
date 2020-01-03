rm(list = ls())

# x <- seq(0, 1, length= 100)
# y <- x
# z <- outer(x, y, prior.fut)
# wireframe(z, drape=T, col.regions=rainbow(100))

require(rmutil)
require(lattice)
require(pracma)

# prior parameters
sigma<-1
sigma.fut<-9.465
alpha<-1
delta.enth<-0.12
delta.skpt<-0
delta.intr<-(delta.skpt+delta.enth)/2
mu<-0.39

## ENTH PRIOR #########################################################
prior.fut<-function(x,y){
  exp(-(sigma*(x-mu)^2)^alpha)*
  exp(-((sigma.fut*((y-x)-delta.skpt))^2)^alpha)
}

 x <- seq(0, 1, length= 100)
 y <- x
 z <- outer(x, y, prior.fut)
 wireframe(z, drape=T, col.regions=rainbow(100))

prior.fut.nc<-tryCatch(integral2(prior.fut,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                        error=function(e) integral2(prior.fut,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])

tail.fut<-tryCatch(
  integral2(prior.fut,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,singular=T)[[1]]/prior.fut.nc,
  error=function(e) 
    integral2(prior.fut,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,singular=T)[[1]]/prior.fut.nc)
tail.fut

tail.fut.intr<-tryCatch(
  integral2(prior.fut,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/prior.fut.nc,
  error=function(e) 
    integral2(prior.fut,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/prior.fut.nc)
tail.fut.intr