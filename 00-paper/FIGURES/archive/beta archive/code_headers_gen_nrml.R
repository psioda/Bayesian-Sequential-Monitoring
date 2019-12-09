# Spike and slab prior information
epsilon<-1E-6
mu0.enth<-(alpha.enth-1)/(alpha.enth+beta.enth-2)

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

mu0.skpt<-(alpha.skpt-1)/(alpha.skpt+beta.skpt-2)

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
