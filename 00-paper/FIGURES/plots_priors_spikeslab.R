## Design parameters, spike/slab case ##
alpha.skpt.1<-0.969595  
beta.skpt.1<-2.186593 
alpha.skpt.2<-10.199189 
beta.skpt.2<-46.344410 
mix.1<-0.1421888
alpha.enth.1<-2.460314 
beta.enth.1<-3.603074 
alpha.enth.2<-96.979737
beta.enth.2<-143.682052
mix.2<-0.3426655 

par(mfrow = c(1, 1)) 
x<-seq(0,1,by=0.01)


plot(x,mix.1*dbeta(x,alpha.skpt.1,beta.skpt.1)+(1-mix.1)*dbeta(x,alpha.skpt.2,beta.skpt.2),
     type="l",xlab="Response Probability",ylab="Density Value",
     main="Spike-Slab Skeptical Prior",
     ylim=c(0,8))
lines(x,dbeta(x,alpha.skpt,beta.skpt))


plot(x,mix.2*dbeta(x,alpha.enth.1,beta.enth.1)+(1-mix.2)*dbeta(x,alpha.enth.2,beta.enth.2),
     type="l",xlab="Response Probability",ylab="Density Value",
     main="Spike-Slab Enthuastic Prior",
     ylim=c(0,10))
lines(x,dbeta(x,alpha.enth,beta.enth))