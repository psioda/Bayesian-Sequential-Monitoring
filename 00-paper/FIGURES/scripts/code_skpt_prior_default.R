mu0.skpt<-p.skpt

## vary sigma to find desired tail area
result<-NA
sigma0.seq<-seq(.01,0.5,by=0.0001)

for (i in 1:length(sigma0.seq)){
  # density & normalizing constant
  sigma0.skpt<-sigma0.seq[i]
  prior.skpt<-function(x){
    exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
  prior.nc.skpt<-function(x){
    exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
  }
  # tail area
  result[i]<-integrate(prior.nc.skpt,lower=p.enth,upper=1-epsilon)[[1]]
}

## find closest match to desired tail area
i<-which(abs(result-tail.skpt)==min(abs(result-tail.skpt)))
sigma0.skpt<-sigma0.seq[i]

# density & normalizing constant
prior.skpt<-function(x){
  exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
}
nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
prior.nc.skpt<-function(x){
  exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
}