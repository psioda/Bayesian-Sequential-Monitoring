## Posterior Mean & Coverage Probability ###################################
pm_cp<-function(index){
  
  y1.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==1)
  y0.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==0)
  y1.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==1)
  y0.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==0)        
  
  posterior.skpt<-function(x,y){
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  
  posterior.skpt.nc<-tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.skpt.nc.x<-function(x,y){
      x*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/
      posterior.skpt.nc
  }
  
  posterior.skpt.nc.y<-function(x,y){
      y*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/
      posterior.skpt.nc
  }
  
  posterior.enth<-function(x,y){
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)
    
  }
  
  posterior.enth.nc<-tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth.nc.x<-function(x,y){
      x*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/posterior.enth.nc
  }
  
  posterior.enth.nc.y<-function(x,y){
      y*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/posterior.enth.nc
  }
  
  mean.enth.x<- tryCatch(tryCatch(integral2(posterior.enth.nc.x,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                error=function(e) integral2(posterior.enth.nc.x,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]]),
                error=function(e) int2(posterior.enth.nc.x,a=c(0,0),b=c(1,1))[[1]])
  mean.enth.y<- tryCatch(tryCatch(integral2(posterior.enth.nc.y,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                error=function(e) integral2(posterior.enth.nc.y,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]]),
                error=function(e) int2(posterior.enth.nc.y,a=c(0,0),b=c(1,1))[[1]])
  mean.skpt.x<- tryCatch(tryCatch(integral2(posterior.skpt.nc.x,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                error=function(e) integral2(posterior.skpt.nc.x,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]]),
                error=function(e) int2(posterior.skpt.nc.x,a=c(0,0),b=c(1,1))[[1]])
  mean.skpt.y<- tryCatch(tryCatch(integral2(posterior.skpt.nc.y,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                error=function(e) integral2(posterior.skpt.nc.y,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]]),
                error=function(e) int2(posterior.skpt.nc.y,a=c(0,0),b=c(1,1))[[1]])
  
  pm.mean.x<-posterior.skpt.nc/(posterior.skpt.nc+posterior.enth.nc)*mean.skpt.x+
             posterior.enth.nc/(posterior.skpt.nc+posterior.enth.nc)*mean.enth.x
  
  pm.mean.y<-posterior.skpt.nc/(posterior.skpt.nc+posterior.enth.nc)*mean.skpt.y+
             posterior.enth.nc/(posterior.skpt.nc+posterior.enth.nc)*mean.enth.y
  
  ## coordinates for test point
  grid.x<-p.PC
  grid.y<-p.IP
  
  ## grid point closest to test point
  grid<-expand.grid(seq(0,1,by=0.01),seq(0,1,by=0.01))
  grid.index<-which.min(sqrt((grid$Var1-grid.x)^2+(grid$Var2-grid.y)^2))
  
  ## probability at cutoff point for credible interval
  final_density<-function(x,y){
      posterior.skpt.nc/(posterior.skpt.nc+posterior.enth.nc)*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/posterior.skpt.nc+
      posterior.enth.nc/(posterior.skpt.nc+posterior.enth.nc)*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/posterior.enth.nc
  }
  
  grid.eval<-final_density(grid$Var1,grid$Var2)/sum(final_density(grid$Var1,grid$Var2))
  grid.eval[is.nan(grid.eval)] <- 0 # 10-29-2019
  grid.index2<-which.min(abs(cumsum(sort(grid.eval))-0.05))
  
  ## is grid point in credible interval?
  coverage<-grid.eval[grid.index]>=sort(grid.eval)[grid.index2]  
  
  return(cbind(pm.mean.x,pm.mean.y,coverage))
}

## Efficacy/Futility Function #############################################
eff_fut<-function(index){
  
  y1.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==1)
  y0.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==0)
  y1.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==1)
  y0.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==0)        
  
  posterior.skpt<-function(x,y){
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  
  posterior.skpt.nc<-tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth<-function(x,y){
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)
  }
  
  posterior.enth.nc<-tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  efficacy<-tryCatch(integral2(posterior.skpt,xmin=0,xmax=1-delta.skpt,ymin=function(x) x+delta.skpt,ymax=1,singular=T)[[1]]/posterior.skpt.nc,
   error=function(e) integral2(posterior.skpt,xmin=0,xmax=1-delta.skpt,ymin=function(x) x+delta.skpt,ymax=1,abstol=1E-6)[[1]]/posterior.skpt.nc)
  
  futility<-tryCatch(1-integral2(posterior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/posterior.enth.nc,
   error=function(e) 1-integral2(posterior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/posterior.enth.nc)
  
  return(cbind(efficacy,futility))
}

## Final Inference ###################################
inference<-function(index){

    y1.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==1)
    y0.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==0)
    y1.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==1)
    y0.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==0)        
    
    posterior.skpt<-function(x,y){
        exp(y1.PC*log(x)+y0.PC*log(1-x))*
        exp(y1.IP*log(y)+y0.IP*log(1-y))*
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
        exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)
    }
    
    posterior.skpt.nc<-tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
              error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
    posterior.enth<-function(x,y){
        exp(y1.PC*log(x)+y0.PC*log(1-x))*
        exp(y1.IP*log(y)+y0.IP*log(1-y))*
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
        exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)
    }
    
    posterior.enth.nc<-tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
              error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
    
    mix.prob<-c(0,0.25,0.5,0.75,1)
    inf.result<-rep(NA,length(mix.prob))
    
    for (m in 1:length(mix.prob)){
    
    omega<-(posterior.skpt.nc*(1-mix.prob[m]))/(posterior.skpt.nc*(1-mix.prob[m])+posterior.enth.nc*mix.prob[m])
    
    final_density<-function(x,y){
      omega*
        exp(y1.PC*log(x)+y0.PC*log(1-x))*
        exp(y1.IP*log(y)+y0.IP*log(1-y))*
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
        exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/posterior.skpt.nc+
      (1-omega)*
        exp(y1.PC*log(x)+y0.PC*log(1-x))*
        exp(y1.IP*log(y)+y0.IP*log(1-y))*
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
        exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/posterior.enth.nc
    }
    
    inf.result[m]<-tryCatch(integral2(final_density,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]],
          error=function(e) integral2(final_density,xmin=0,xmax=1,ymin=function(x) x,ymax=1,abstol=1E-6)[[1]])
    }
    
    return(inf.result>sig.eff)
}

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
  
  assign("mu",mu,envir = .GlobalEnv)
  assign("sigma0.placebo",sigma0,envir = .GlobalEnv)
  assign("lambda0.placebo",lambda0,envir = .GlobalEnv)
  
  prior.placebo<-function(x){exp(-(abs(x-mu)/sigma0)^lambda0)}
  nc<-integrate(prior.placebo,lower=0+epsilon,upper=1-epsilon)[[1]]
  prior.placebo.nc<-function(x){exp(-(abs(x-mu)/sigma0)^lambda0)/nc}
  
  return(prior.placebo.nc)
}

skpt_tail_area<-function(){
  
  sigma0.seq<-seq(.01,0.5,by=0.01)
  lambda0.seq<-seq(1,2,by=0.1)
  result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  q.outer<-0.025  # y > x + delta.enth
  q.inner<-0.010  # y > x + delta.intr
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.skpt<-sigma0.seq[i]
      lambda0.skpt<-lambda0.seq[j]
      
      prior.skpt<-function(x,y){
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
          exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)}
      nc<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
 error=function(e) integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
      
      result1[i,j]<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,singular=T)[[1]]/nc,
                             error=function(e) integral2(prior.skpt,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,abstol=1E-6)[[1]]/nc)
      result2[i,j]<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/nc, 
                             error=function(e) integral2(prior.skpt,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/nc)
    }
  }
  
  result3=abs(result1-q.outer)+abs(result2-q.inner)
  index<-which(result3 == min(result3), arr.ind = TRUE)
  i<-index[1]
  j<-index[2]
  sigma0.skpt<-sigma0.seq[i]
  lambda0.skpt<-lambda0.seq[j]
  print(result1[i,j])
  print(result2[i,j])
  
  assign("sigma0.skpt",sigma0.skpt,envir = .GlobalEnv)
  assign("lambda0.skpt",lambda0.skpt,envir = .GlobalEnv)
  
  prior.skpt<-function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)}
  nc<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
               error=function(e) integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  prior.skpt.nc<-function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/nc}
  
}

enth_tail_area<-function(){
  
  sigma0.seq<-seq(0.085,0.90,by=0.001)
  lambda0.seq<-seq(2,2,by=0.1)
  result1<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  #result2<-matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  q.outer<-0.025  # y > x + delta.enth
  #q.inner<-0.010  # y > x + delta.intr
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.enth<-sigma0.seq[i]
      lambda0.enth<-lambda0.seq[j]
      
      prior.enth<-function(x,y){
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
          exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)}
      nc<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
 error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
      
      result1[i,j]<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]]/nc,
                             error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,abstol=1E-6)[[1]]/nc)      
      #result2[i,j]<-tryCatch(integral2(prior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/prior.enth.nc,
      #     error=function(e) integral2(prior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/prior.enth.nc)
    }
  }
  # 1-result1 since looking for q.outer area that y < x
  # 1-result2 since looking for q.inner area that y < x + delta.intr
  result3=abs(1-result1-q.outer)#+abs(1-result2-q.inner)
  index<-which(result3 == min(result3), arr.ind = TRUE)
  i<-index[1]
  j<-index[2]
  sigma0.enth<-sigma0.seq[i]
  lambda0.enth<-lambda0.seq[j]
  # 1-result1 since looking for q.outer area that y < x
  # 1-result2 since looking for q.inner area that y < x + delta.intr
  print(1-result1[i,j]-q.outer)
  #print(1-result2[i,j])
  
  assign("sigma0.enth",sigma0.enth,envir = .GlobalEnv)
  assign("lambda0.enth",lambda0.enth,envir = .GlobalEnv)
  
  prior.enth<-function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)}
  nc<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  prior.enth.nc<-function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/nc}
  
  
  return(prior.enth.nc)
}