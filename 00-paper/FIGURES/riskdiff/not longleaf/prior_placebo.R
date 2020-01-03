fcn_prior_placebo<-function(){
  
sigma0.seq<-seq(.15,35,by=0.01)
lambda0.seq<-seq(4.5,5.5,by=0.1)
result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))

q.outer<-0.5  # window around mu +/- delta.enth
q.inner<-0.25 # window around mu +/- delta.intr

epsilon<-0

  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0<-sigma0.seq[i]
      lambda0<-lambda0.seq[j]
      
      prior.placebo<-function(x){exp(-(abs(x-mu)/sigma0)^lambda0)}
      nc<-integrate(prior.placebo,lower=0+epsilon,upper=1-epsilon)[[1]]
      prior.placebo.nc<-function(x){exp(-(abs(x-mu)/sigma0)^lambda0)/nc}
      
      result1[i,j]<-integrate(prior.placebo.nc,lower=mu-delta.enth,upper=mu+delta.enth)[[1]]
      result2[i,j]<-integrate(prior.placebo.nc,lower=mu-delta.intr,upper=mu+delta.intr)[[1]]
    }
  }

  result3=abs(result1-q.outer)+abs(result2-q.inner)
  index<-which(result3 == min(result3), arr.ind = TRUE)
  i<-index[1]
  j<-index[2]
  print(result1[i,j])
  print(result2[i,j])
  sigma0<-sigma0.seq[i]
  lambda0<-lambda0.seq[j]
  
  xx<-seq(0,1,length=1000)  
  plot(xx,prior.placebo.nc(xx),ylim=c(0,max(prior.placebo.nc(xx))),type='l')
  
  assign("mu",mu,envir = .GlobalEnv)
  assign("sigma0.placebo",sigma0,envir = .GlobalEnv)
  assign("lambda0.placebo",lambda0,envir = .GlobalEnv)
}