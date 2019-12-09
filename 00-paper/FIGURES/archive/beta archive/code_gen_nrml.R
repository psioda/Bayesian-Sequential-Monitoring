rm(list = ls())

## 3 parameters for generalized normal distribution
mu0<-0.4
sigma0<-1
lambda0<-2.5 # does not change height

## vary sigma to find desired tail area
result<-NA
sigma0.seq<-seq(.01,0.5,by=0.0001)
epsilon<-1E-6
for (i in 1:length(sigma0.seq)){
# density & normalizing constant
sigma0<-sigma0.seq[i]
prior<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)
}
nc<-integrate(prior,lower=0+epsilon,upper=1-epsilon)[[1]]
prior.nc<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)/nc
}
# tail area
result[i]<-integrate(prior.nc,lower=epsilon,upper=0.2)[[1]]
}

## find closest match to desired tail area
i<-which(abs(result-0.05)==min(abs(result-0.05)))
sigma0<-sigma0.seq[i]

# density & normalizing constant
prior<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)
}
nc<-integrate(prior,lower=0+epsilon,upper=1-epsilon)[[1]]
prior.nc<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)/nc
}

# plot prior
x<-seq(0.01,0.99,by=0.01)
plot(x,prior.nc(x),type='l')
max(prior.nc(x))
integrate(prior.nc,lower=epsilon,upper=0.2)[[1]]

# expected value
prior.nc.exp<-function(x){
  x*exp(-(abs(x-mu0)/sigma0)^lambda0)/nc
}
integrate(prior.nc.exp,lower=epsilon,upper=1-epsilon)[[1]]
# add data
n<-5
y<-2

posterior<-function(x){
  exp(y*log(x)+(n-y)*log(1-x)-((abs(x-mu0)/sigma0)^lambda0))
}
nc<-integrate(posterior,lower=0+epsilon,upper=1-epsilon)[[1]]
posterior.nc<-function(x){
exp(y*log(x)+(n-y)*log(1-x)-((abs(x-mu0)/sigma0)^lambda0))/nc
}

integrate(posterior.nc,lower=0.2,upper=1-epsilon)[[1]]

