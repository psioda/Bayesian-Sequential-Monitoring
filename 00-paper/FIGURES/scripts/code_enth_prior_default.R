mu0.enth<-p.enth

## vary sigma to find desired tail area
result<-NA
sigma0.seq<-seq(.01,0.5,by=0.0001)

for (i in 1:length(sigma0.seq)){
  # density & normalizing constant
  sigma0.enth<-sigma0.seq[i]
  prior.enth<-function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
  }
  nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
  prior.nc.enth<-function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
  }
  # tail area
  result[i]<-integrate(prior.nc.enth,lower=epsilon,upper=p.skpt)[[1]]
}

## find closest match to desired tail area
i<-which(abs(result-tail.enth)==min(abs(result-tail.enth)))
sigma0.enth<-sigma0.seq[i]

# density & normalizing constant
prior.enth<-function(x){
  exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
}
nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
prior.nc.enth<-function(x){
  exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
}