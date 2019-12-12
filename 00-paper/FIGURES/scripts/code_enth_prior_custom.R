# Spike and slab prior information
mu0.enth<-p.enth

## vary sigma to find desired tail area
sigma0.seq<-seq(.01,2,by=0.01)
lambda0.seq<-seq(0.1,7,by=0.1)
result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))

for (i in 1:length(sigma0.seq)){
  for (j in 1:length(lambda0.seq)){
    # density & normalizing constant
    sigma0.enth<-sigma0.seq[i]
    lambda0.enth<-lambda0.seq[j]
    prior.enth<-function(x){
      exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
    }
    nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
    prior.nc.enth<-function(x){
      exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
    }
    # tail area
    result1[i,j]<-integrate(prior.nc.enth,lower=0+epsilon,upper=p.skpt)[[1]]
    result2[i,j]<-integrate(prior.nc.enth,lower=p.intr,upper=p.enth)[[1]]
  }
}
result3=abs(result1-tail.enth)+abs(result2-0.3396372*scale)
index<-which(result3 == min(result3), arr.ind = TRUE)

## Check work
i<-index[1]
j<-index[2]
sigma0.enth<-sigma0.seq[i]
lambda0.enth<-lambda0.seq[j]
prior.enth<-function(x){
  exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
}
nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
prior.nc.enth<-function(x){
  exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
}
sigma0.enth
lambda0.enth
integrate(prior.nc.enth,lower=0+epsilon,upper=p.skpt)[[1]]
integrate(prior.nc.enth,lower=p.intr,upper=p.enth)[[1]]
