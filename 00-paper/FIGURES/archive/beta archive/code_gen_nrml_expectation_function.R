y1[n[l]]<-20
y0[n[l]]<-20
######################
### SKEPTICAL CASE ###
######################
prior.skpt<-function(x){
  exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
}
prior.skpt.nc<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]

posterior.skpt<-function(x){
  exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
}
posterior.skpt.nc<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]

post.prob.skpt<-posterior.skpt.nc/prior.skpt.nc

posterior.skpt.nc.exp<-function(x){
  x*exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/posterior.skpt.nc
}
skpt.exp<-integrate(posterior.skpt.nc.exp,lower=0+epsilon,upper=1-epsilon)[[1]]

#######################
### ENTHUASTIC CASE ###
#######################
prior.enth<-function(x){
  exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
}
prior.enth.nc<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]

posterior.enth<-function(x){
  exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
}
posterior.enth.nc<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]

post.prob.enth<-posterior.enth.nc/prior.enth.nc

posterior.enth.nc.exp<-function(x){
  x*exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/posterior.enth.nc
}
enth.exp<-integrate(posterior.enth.nc.exp,lower=0+epsilon,upper=1-epsilon)[[1]]

####################
### WEIGHTED AVG ###
####################
c<-post.prob.skpt/(post.prob.skpt+post.prob.enth)
print(c*skpt.exp+(1-c)*enth.exp)

