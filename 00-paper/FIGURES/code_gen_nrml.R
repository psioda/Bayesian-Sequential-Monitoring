rm(list = ls())

mu0<-0.4
sigma0<-1
lambda0<-4 # does not change height
epsilon<-1E-6

result<-NA
sigma0.seq<-seq(.1,0.5,by=0.01)

for (i in 1:length(sigma0.seq)){
sigma0<-sigma0.seq[i]
  
gen.nrml<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)
}
nc<-integrate(gen.nrml,lower=0+epsilon,upper=1-epsilon)[[1]]
gen.nrml.nc<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)/nc
}

x<-seq(0.01,0.99,by=0.01)
plot(x,gen.nrml.nc(x),type='l')
max(gen.nrml.nc(x))
#integrate(gen.nrml.nc,lower=epsilon,upper=1-epsilon)
result[i]<-integrate(gen.nrml.nc,lower=epsilon,upper=0.2)[[1]]
}

min(abs(result-0.05))
which(result==c(min(abs(result-0.05))+0.05))

i<-22
sigma0<-sigma0.seq[i]

gen.nrml<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)
}
nc<-integrate(gen.nrml,lower=0+epsilon,upper=1-epsilon)[[1]]
gen.nrml.nc<-function(x){
  exp(-(abs(x-mu0)/sigma0)^lambda0)/nc
}

x<-seq(0.01,0.99,by=0.01)
plot(x,gen.nrml.nc(x),type='l')
max(gen.nrml.nc(x))
#integrate(gen.nrml.nc,lower=epsilon,upper=1-epsilon)
result[i]<-integrate(gen.nrml.nc,lower=epsilon,upper=0.2)[[1]]
