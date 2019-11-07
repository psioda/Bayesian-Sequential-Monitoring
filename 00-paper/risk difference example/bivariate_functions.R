## Posterior Mean & Coverage Probability ###################################
pm_cp<-function(index){
  
  # need to get model probabilities from skeptical and enthuastic prior
  # need to get expected value from skeptical and enthuastic prior
  # need to weight them
  
  y1.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==1)
  y0.IP<-sum(responses.IP[outcome.times.IP<=outcome.times.all[index]]==0)
  y1.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==1)
  y0.PC<-sum(responses.PC[outcome.times.PC<=outcome.times.all[index]]==0)        
  
  posterior.skpt<-function(x,y){
    exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.skpt*(x-mu)^2)^alpha)*
      exp(-((sigma.skpt*((y-x)-delta.skpt))^2)^alpha)
  }
  
  
  posterior.skpt.nc<-tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                              error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.skpt.nc.x<-function(x,y){
    x*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.skpt*(x-mu)^2)^alpha)*
      exp(-((sigma.skpt*((y-x)-delta.skpt))^2)^alpha)/
      posterior.skpt.nc
  }
  
  posterior.skpt.nc.y<-function(x,y){
    y*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.skpt*(x-mu)^2)^alpha)*
      exp(-((sigma.skpt*((y-x)-delta.skpt))^2)^alpha)/
      posterior.skpt.nc
  }
  
  posterior.enth<-function(x,y){
    exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.enth*(x-mu)^2)^alpha)*
      exp(-((sigma.enth*((y-x)-delta.enth))^2)^alpha)
  }
  
  posterior.enth.nc<-tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                              error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth.nc.x<-function(x,y){
    x*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.enth*(x-mu)^2)^alpha)*
      exp(-((sigma.enth*((y-x)-delta.enth))^2)^alpha)/
      posterior.enth.nc
  }
  
  posterior.enth.nc.y<-function(x,y){
    y*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.enth*(x-mu)^2)^alpha)*
      exp(-((sigma.enth*((y-x)-delta.enth))^2)^alpha)/
      posterior.enth.nc
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
      exp(-(sigma.skpt*(x-mu)^2)^alpha)*
      exp(-((sigma.skpt*((y-x)-delta.skpt))^2)^alpha)/posterior.skpt.nc+
      posterior.enth.nc/(posterior.skpt.nc+posterior.enth.nc)*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(sigma.enth*(x-mu)^2)^alpha)*
      exp(-((sigma.enth*((y-x)-delta.enth))^2)^alpha)/posterior.enth.nc
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
      exp(-(sigma.skpt*(x-mu)^2)^alpha)*
      exp(-((sigma.skpt*((y-x)-delta.skpt))^2)^alpha)
  }
  
  posterior.skpt.nc<-tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                              error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth<-function(x,y){
    exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(-(sigma.enth*(x-mu)^2)^alpha)*
      exp(-((sigma.enth*((y-x)-delta.enth))^2)^alpha)
  }
  posterior.enth.nc<-tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                              error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  efficacy<-tryCatch(
    integral2(posterior.skpt,xmin=0,xmax=1-delta.skpt,
              ymin=function(x) x+delta.skpt,ymax=1,singular=T)[[1]]/posterior.skpt.nc,
    error=function(e) 
      integral2(posterior.skpt,xmin=0,xmax=1-delta.skpt,
                ymin=function(x) x+delta.skpt,ymax=1,abstol=1E-6)[[1]]/posterior.skpt.nc)
  
  futility<-tryCatch(
    1-integral2(posterior.enth,xmin=0,xmax=1-delta.intr,
                ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/posterior.enth.nc,
    error=function(e) 
      1-integral2(posterior.enth,xmin=0,xmax=1-delta.intr,
                  ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/posterior.enth.nc)
  
  
  return(cbind(efficacy,futility))
}