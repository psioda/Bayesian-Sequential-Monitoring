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
  futility<-integrate(posterior.nc.enth,lower=0+epsilon,upper=p.intr)[[1]]
  
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
  
  return(cbind(pm.mean,inf.result>sig.eff,p.range[j]>=lower_cr & p.range[j]<=upper_cr,lower_cr,upper_cr))
}