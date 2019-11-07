rm(list = ls())

# x <- seq(0, 1, length= 100)
# y <- x
# z <- outer(x, y, prior.fut)
# wireframe(z, drape=T, col.regions=rainbow(100))

require(rmutil)
require(lattice)
require(pracma)

# prior parameters
sigma<-9.435
alpha<-2
delta.skpt<-0
mu<-0.1

## ENTH PRIOR #########################################################
prior.fut<-function(x,y){
  exp(-(sigma*(x-mu)^2)^alpha)*
  exp(-((sigma*((y-x)-delta.skpt))^2)^alpha)
}

 x <- seq(0, 1, length= 100)
 y <- x
 z <- outer(x, y, prior.fut)
 wireframe(z, drape=T, col.regions=rainbow(100))

prior.fut.nc<-tryCatch(integral2(prior.fut,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                        error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])

tail.fut<-tryCatch(
  integral2(prior.fut,xmin=0,xmax=0.9,ymin=function(x) x+0.1,ymax=1,singular=T)[[1]]/prior.fut.nc,
  error=function(e) 
    integral2(prior.fut,xmin=0,xmax=0.9,ymin=function(x) x+0.1,ymax=1,singular=T)[[1]]/prior.fut.nc)
tail.fut