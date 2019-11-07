# x <- seq(0, 1, length= 100)
# y <- x
# z <- outer(x, y, prior.enth)
# wireframe(z, drape=T, col.regions=rainbow(100))

require(rmutil)
require(lattice)
require(pracma)

# prior parameters
sigma<-9.5
alpha<-2
delta.enth<-0.1
delta.skpt<-0
delta.intr=(delta.skpt+delta.enth)/2
mu<-0.2

## ENTH PRIOR #########################################################
prior.enth<-function(x,y){
  #exp(-(sigma*(x-mu)^2)^alpha)*
  exp(-((sigma*((y-x)-delta.enth))^2)^alpha)
}

# x <- seq(0, 1, length= 100)
# y <- x
# z <- outer(x, y, prior.enth)
# wireframe(z, drape=T, col.regions=rainbow(100))

prior.enth.nc<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                        error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])

tail.enth<-tryCatch(
  integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]]/prior.enth.nc,
  error=function(e) 
    integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]]/prior.enth.nc)
1-tail.enth