rm(list = ls())
require(rmutil)
require(lattice)
require(pracma)

sigma<-4
alpha<-2
delta.enth<-0.2
delta.skpt<-0

# sample data: stop early for futility
y1.x<-seq(0,20,by=2)
y0.x<-seq(0,30,by=3)
y1.y<-seq(0,20,by=2)
y0.y<-seq(0,30,by=3)

# # sample data: stop early for efficacy
# y1.x<-c(0,1,2,3,4,5,6,7)
# y0.x<-c(0,4,8,12,16,20,24,28)
# y1.y<-c(0,3,6,9,12,15,18,21)
# y0.y<-c(0,2,4,6,8,10,12,14)

eff<-NA
fut<-NA

for (i in 1:length(y1.x)){
  
  posterior.skpt<-function(x,y){
    exp(y1.x[i]*log(x)+y0.x[i]*log(1-x))*
      exp(y1.y[i]*log(y)+y0.y[i]*log(1-y))*
      exp(-((sigma*((y-x)-delta.skpt))^2)^alpha)
  }
  posterior.skpt.nc<-integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1)[[1]]
  
  posterior.enth<-function(x,y){
    exp(y1.x[i]*log(x)+y0.x[i]*log(1-x))*
      exp(y1.y[i]*log(y)+y0.y[i]*log(1-y))*
      exp(-((sigma*((y-x)-delta.enth))^2)^alpha)
  }
  posterior.enth.nc<-integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1)[[1]]

  ymax <- function(x) x+0.1
  eff[i]<-1-integral2(posterior.skpt,xmin=0,xmax=.9,ymin=.1,ymax)[[1]]/posterior.skpt.nc
  
  ymax <- function(x) x+0.15
  fut[i]<-integral2(posterior.enth,xmin=0,xmax=.85,ymin=.15,ymax)[[1]]/posterior.enth.nc
}
eff
fut