# Spike and slab prior information
mu0.skpt<-p.skpt

## vary sigma to find desired tail area
sigma0.seq<-seq(.01,2,by=0.01)
lambda0.seq<-seq(0.1,7,by=0.1)
result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))

for (i in 1:length(sigma0.seq)){
  for (j in 1:length(lambda0.seq)){
    # density & normalizing constant
    sigma0.skpt<-sigma0.seq[i]
    lambda0.skpt<-lambda0.seq[j]
    prior.skpt<-function(x){
      exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
    }
    nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
    prior.nc.skpt<-function(x){
      exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
    }
    # tail area
    result1[i,j]<-integrate(prior.nc.skpt,lower=p.enth,upper=1-epsilon)[[1]]
    result2[i,j]<-integrate(prior.nc.skpt,lower=p.skpt,upper=p.intr)[[1]]
  }
}
result3=abs(result1-tail.skpt)+abs(result2-0.3370991*scale)
index<-which(result3 == min(result3), arr.ind = TRUE)

## Check work
i<-index[1]
j<-index[2]
sigma0.skpt<-sigma0.seq[i]
lambda0.skpt<-lambda0.seq[j]
prior.skpt<-function(x){
  exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
}
nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
prior.nc.skpt<-function(x){
  exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
}
sigma0.skpt
lambda0.skpt
integrate(prior.nc.skpt,p.enth,1)
integrate(prior.nc.skpt,p.skpt,p.intr)
