## Efficacy/Futility Function #############################################
eff_fut<-function(index){
  
  posterior.skpt<-function(x){
    exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
  }
  skpt.nc<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
  posterior.nc.skpt<-function(x){
    exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/skpt.nc
  }
  efficacy<-integrate(posterior.nc.skpt,lower=0+epsilon,upper=p.skpt)[[1]]
  
  posterior.enth<-function(x){
    exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
  }
  enth.nc<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
  posterior.nc.enth<-function(x){
    exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/enth.nc
  }
  futility<-integrate(posterior.nc.enth,lower=0+epsilon,upper=p.enth)[[1]] # changed on 2020-09-09 upper = p.intr changed to upper = p.enth
  
  return(cbind(efficacy,futility))
}

## Posterior Mean, Coverage Probability, Final Inference ##################
pm_cp<-function(index){
  
  posterior.skpt<-function(x){
    exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
  }
  skpt.nc<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
 
  posterior.enth<-function(x){
    exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
  }
  enth.nc<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
  
  omega<-skpt.nc/(skpt.nc+enth.nc)
  
  final_density<-function(x){
        omega*exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/skpt.nc+
    (1-omega)*exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/enth.nc
  }
  
  final_density_x<-function(x){
    x*(
          omega*exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/skpt.nc+
      (1-omega)*exp(y1[index]*log(x)+y0[index]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/enth.nc
      )
  }
  
  inf.result<-integrate(final_density,lower=p.skpt,upper=1-epsilon)[[1]]
  pm.mean<-integrate(final_density_x,lower=0+epsilon,upper=1-epsilon)[[1]]
 
  lower_cr<-0+1e-3
  tail <-0
  while(tail<cred.tail/2){
    tail<-integrate(final_density,0,lower_cr)[[1]]
    if (tail<=cred.tail/2) lower_cr<-lower_cr+1e-3
  }
  
  upper_cr<-1-1e-3 # added 12/8/19
  tail <-0
  while(tail<cred.tail/2){
    tail<-integrate(final_density,upper_cr,1)[[1]]
    if (tail<=cred.tail/2) upper_cr<-upper_cr-1e-3
  }
  
  return(cbind(pm.mean,
               inf.result>sig.eff,
               tryCatch(p.range[j]>=lower_cr & p.range[j]<=upper_cr,error=function(e) NA),
               lower_cr,
               upper_cr))
}

skpt_prior_default<-function(){
  
  mu0.skpt<-p.skpt
  sigma0.seq<-seq(.01,0.5,by=0.0001)
  lambda0.skpt<-2
  
  result<-NA
  
  for (i in 1:length(sigma0.seq)){
    
    sigma0.skpt<-sigma0.seq[i]
    
    prior.skpt<-function(x){
      exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
    }
    nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
    prior.nc.skpt<-function(x){
      exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
    }
    
    result[i]<-integrate(prior.nc.skpt,lower=p.enth,upper=1-epsilon)[[1]]
  }
  
  i<-which(abs(result-tail.skpt)==min(abs(result-tail.skpt)))
  sigma0.skpt<-sigma0.seq[i]
  
  prior.skpt<-function(x){
    exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
  prior.nc.skpt<-function(x){
    exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
  }
  
  print(paste0("mu: ",mu0.skpt,", sigma: ",sigma0.skpt,", lambda: ",lambda0.skpt))
  print(paste0("Tail area: ",result[i]))
  print(paste0("Half-width area: ",
               integrate(prior.nc.skpt,lower=p.skpt,upper=p.intr)[[1]]))
  
  assign("mu0.skpt",mu0.skpt,envir = .GlobalEnv)
  assign("sigma0.skpt",sigma0.skpt,envir = .GlobalEnv)
  assign("lambda0.skpt",lambda0.skpt,envir = .GlobalEnv)
  
  return(prior.nc.skpt)
}

enth_prior_default<-function(){
  
  mu0.enth<-p.enth
  sigma0.seq<-seq(.01,0.5,by=0.0001)
  lambda0.enth<-2
  
  result<-NA
  
  for (i in 1:length(sigma0.seq)){
    
    sigma0.enth<-sigma0.seq[i]
    
    prior.enth<-function(x){
      exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
    }
    nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
    prior.nc.enth<-function(x){
      exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
    }
    
    result[i]<-integrate(prior.nc.enth,lower=0+epsilon,upper=p.skpt)[[1]]
  }
  
  i<-which(abs(result-tail.enth)==min(abs(result-tail.enth)))
  sigma0.enth<-sigma0.seq[i]
  
  prior.enth<-function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
  }
  nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
  prior.nc.enth<-function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
  }
  
  print(paste0("mu: ",mu0.enth,", sigma: ",sigma0.enth,", lambda: ",lambda0.enth))
  print(paste0("Tail area: ",result[i]))
  print(paste0("Half-width area: ",
               integrate(prior.nc.enth,lower=p.intr,upper=p.enth)[[1]]))
  
  assign("mu0.enth",mu0.enth,envir = .GlobalEnv)
  assign("sigma0.enth",sigma0.enth,envir = .GlobalEnv)
  assign("lambda0.enth",lambda0.enth,envir = .GlobalEnv)
  
  return(prior.nc.enth)
}

skpt_prior_custom<-function(scale){
  
  mu0.skpt<-p.skpt
  sigma0.seq<-seq(.01,2,by=0.001)
  lambda0.seq<-seq(0.01,2,by=0.01)
  result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.skpt<-sigma0.seq[i]
      lambda0.skpt<-lambda0.seq[j]
      
      prior.skpt<-function(x){
        exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
      }
      nc.skpt<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
      prior.nc.skpt<-function(x){
        exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt
      }
      
      result1[i,j]<-integrate(prior.nc.skpt,lower=p.enth,upper=1-epsilon)[[1]]
      result2[i,j]<-integrate(prior.nc.skpt,lower=p.intr,upper=p.enth)[[1]]
    }
  }
  result3=abs(result1-tail.skpt)+abs(result2-(pnorm(qnorm(tail.skpt)/2) - tail.skpt)*scale)
  index<-which(result3 == min(result3), arr.ind = TRUE)
  
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
  
  print(paste0("mu: ",mu0.skpt,", sigma: ",sigma0.skpt,", lambda: ",lambda0.skpt))
  print(paste0("Tail area: ",result1[i,j]))
  print(paste0("Half-width area: ",result2[i,j]))
  
  assign("mu0.skpt",mu0.skpt,envir = .GlobalEnv)
  assign("sigma0.skpt",sigma0.skpt,envir = .GlobalEnv)
  assign("lambda0.skpt",lambda0.skpt,envir = .GlobalEnv)
  
  return(prior.nc.skpt)
}

enth_prior_custom<-function(scale){
  
  mu0.enth<-p.enth
  sigma0.seq<-seq(.01,0.5,by=0.001)
  lambda0.seq<-seq(2,7,by=0.01)
  result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.enth<-sigma0.seq[i]
      lambda0.enth<-lambda0.seq[j]
      
      prior.enth<-function(x){
        exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
      }
      nc.enth<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
      prior.nc.enth<-function(x){
        exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
      }
      
      result1[i,j]<-integrate(prior.nc.enth,lower=0+epsilon,upper=p.skpt)[[1]]
      result2[i,j]<-integrate(prior.nc.enth,lower=p.skpt,upper=p.intr)[[1]]
    }
  }
  result3=abs(result1-tail.enth)+abs(result2 - (pnorm(qnorm(tail.enth)/2) - tail.enth)*scale)
  index<-which(result3 == min(result3), arr.ind = TRUE)
  
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
  
  print(paste0("mu: ",mu0.enth,", sigma: ",sigma0.enth,", lambda: ",lambda0.enth))
  print(paste0("Tail area: ",result1[i,j]))
  print(paste0("Half-width area: ",result2[i,j]))
  
  assign("mu0.enth",mu0.enth,envir = .GlobalEnv)
  assign("sigma0.enth",sigma0.enth,envir = .GlobalEnv)
  assign("lambda0.enth",lambda0.enth,envir = .GlobalEnv)
  assign("tail.enth.actual",result1[i,j],envir = .GlobalEnv)
  
  return(prior.nc.enth)
}