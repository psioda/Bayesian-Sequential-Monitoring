par(ask=FALSE)
par(mfrow=c(2,2))
#par(mfrow = c(1, 2)) 
x<-seq(0,1,by=0.01)
plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",
     xlab="Response Probability",ylab="Density Value",
     ylim=c(0,7.5))
## Design parameters, spike/slab case
alpha.skpt.1<-0.969595
beta.skpt.1<-2.186593
alpha.skpt.2<-10.199189
beta.skpt.2<-46.344410
mix.1<-0.1421888
lines(x,mix.1*dbeta(x,alpha.skpt.1,beta.skpt.1)+(1-mix.1)*dbeta(x,alpha.skpt.2,beta.skpt.2),
      lty='longdash')

plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",
     xlab="Response Probability",ylab="Density Value",
     ylim=c(0,4.5))
## Design parameters, heavy case
alpha.skpt.1<-4.360271 
beta.skpt.1<-17.626647
alpha.skpt.2<-1.242130 
beta.skpt.2<-4.901135
mix.1<-0.4774422 
lines(x,mix.1*dbeta(x,alpha.skpt.1,beta.skpt.1)+(1-mix.1)*dbeta(x,alpha.skpt.2,beta.skpt.2),
      lty='longdash')

plot(x,dbeta(x,alpha.enth,beta.enth),type="l",
     xlab="Response Probability",ylab="Density Value",
     ylim=c(0,9))
## Design parameters, spike/slab case
alpha.enth.1<-2.460314
beta.enth.1<-3.603074
alpha.enth.2<-96.979737
beta.enth.2<-143.682052
mix.2<-0.3426655
lines(x,mix.2*dbeta(x,alpha.enth.1,beta.enth.1)+(1-mix.2)*dbeta(x,alpha.enth.2,beta.enth.2),
      lty='longdash')

plot(x,dbeta(x,alpha.enth,beta.enth),type="l",
     xlab="Response Probability",ylab="Density Value",
     ylim=c(0,4))
## Design parameters, spike/slab case
alpha.enth.1<-3.810989 
beta.enth.1<-4.908287
alpha.enth.2<-36.729606 
beta.enth.2<-97.651133 
mix.2<-0.7709653 
lines(x,mix.2*dbeta(x,alpha.enth.1,beta.enth.1)+(1-mix.2)*dbeta(x,alpha.enth.2,beta.enth.2),
      lty='longdash')
