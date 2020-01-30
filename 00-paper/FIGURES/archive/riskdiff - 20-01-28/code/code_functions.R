## Posterior Mean & Coverage Probability ###################################
pm_cp <- function(index){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)        
  
  posterior.skpt <- function(x,y){
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  posterior.skpt.nc <- tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth <- function(x,y){
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)
    
  }
  posterior.enth.nc <- tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  omega <- (posterior.skpt.nc*(1-mix.prob))/(posterior.skpt.nc*(1-mix.prob)+posterior.enth.nc*mix.prob)
  final_density <- function(x,y){
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
  
  final_density_x <- function(x,y){
    x*(
      omega*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/posterior.skpt.nc+
      (1-omega)*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/posterior.enth.nc)
  }
  final_density_y <- function(x,y){
    y*(
      omega*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/posterior.skpt.nc+
      (1-omega)*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/posterior.enth.nc)
  }
  
  pm.mean.x <-  tryCatch(tryCatch(integral2(final_density_x,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
              error=function(e) integral2(final_density_x,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]]),
                   error=function(e) int2(final_density_x,a=c(0,0),b=c(1,1))[[1]])
  pm.mean.y <-  tryCatch(tryCatch(integral2(final_density_y,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
              error=function(e) integral2(final_density_y,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]]),
                   error=function(e) int2(final_density_y,a=c(0,0),b=c(1,1))[[1]])
  
  # grid.index returns a 1d index for the test point
  grid.x <- p.PC
  grid.y <- p.IP
  grid <- expand.grid(seq(0,1,by=0.01),seq(0,1,by=0.01))
  grid.index <- which.min(sqrt((grid$Var1-grid.x)^2+(grid$Var2-grid.y)^2))
  # evaluate "normalized density" at the grid
  grid.eval <- final_density(grid$Var1,grid$Var2)/sum(final_density(grid$Var1,grid$Var2))
  grid.eval[is.nan(grid.eval)]  <-  0 # 10-29-2019
  # find the cutoff for cred tail percentile (absolute value unnecessary)
  grid.index2 <- which.min(abs(cumsum(sort(grid.eval))-cred.tail))
  # is grid point in credible interval?
  coverage <- grid.eval[grid.index]>=sort(grid.eval)[grid.index2]  
  
  return(cbind(pm.mean.x,pm.mean.y,coverage))
}

## Final Inference ###################################
inference.skpt <- function(index){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)        
  
  posterior.skpt <- function(x,y){
    exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  posterior.skpt.nc <- tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                              error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth <- function(x,y){
    exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)
  }
  posterior.enth.nc <- tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
                              error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  omega <- (posterior.skpt.nc*(1-mix.prob))/(posterior.skpt.nc*(1-mix.prob)+posterior.enth.nc*mix.prob)
  final_density <- function(x,y){
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
  
  result <- tryCatch(integral2(final_density,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]],
                   error=function(e) integral2(final_density,xmin=0,xmax=1,ymin=function(x) x,ymax=1,abstol=1E-6)[[1]])
  
  return(result)
}

## Efficacy/Futility Function #############################################
eff_fut <- function(index){
  
  y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)        
  
  posterior.skpt <- function(x,y){
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)
  }
  posterior.skpt.nc <- tryCatch(integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  posterior.enth <- function(x,y){
      exp(y1.IP*log(y)+y0.IP*log(1-y))*
      exp(y1.PC*log(x)+y0.PC*log(1-x))*
      exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)
  }
  posterior.enth.nc <- tryCatch(integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(posterior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  
  efficacy <- tryCatch(integral2(posterior.skpt,xmin=0,xmax=1-delta.skpt,ymin=function(x) x+delta.skpt,ymax=1,singular=T)[[1]]/posterior.skpt.nc,
   error=function(e) integral2(posterior.skpt,xmin=0,xmax=1-delta.skpt,ymin=function(x) x+delta.skpt,ymax=1,abstol=1E-6)[[1]]/posterior.skpt.nc)
  
  futility <- tryCatch(1-integral2(posterior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/posterior.enth.nc,
   error=function(e) 1-integral2(posterior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/posterior.enth.nc)
  
  return(cbind(efficacy,futility))
}

fcn_prior_placebo <- function(){
  
  sigma0.seq <- seq(.15,35,by=0.01)
  lambda0.seq <- seq(4.5,5.5,by=0.1)
  result1 <- matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  result2 <- matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  q.outer <- 0.5  # window around mu +/- delta.enth
  q.inner <- 0.25 # window around mu +/- delta.intr
  
  epsilon <- 0
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0 <- sigma0.seq[i]
      lambda0 <- lambda0.seq[j]
      
      prior.placebo <- function(x){exp(-(abs(x-mu)/sigma0)^lambda0)}
      nc <- integrate(prior.placebo,lower=0+epsilon,upper=1-epsilon)[[1]]
      prior.placebo.nc <- function(x){exp(-(abs(x-mu)/sigma0)^lambda0)/nc}
      
      result1[i,j] <- integrate(prior.placebo.nc,lower=mu-delta.enth,upper=mu+delta.enth)[[1]]
      result2[i,j] <- integrate(prior.placebo.nc,lower=mu-delta.intr,upper=mu+delta.intr)[[1]]
    }
  }
  
  result3=abs(result1-q.outer)+abs(result2-q.inner)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  print(result1[i,j])
  print(result2[i,j])
  sigma0 <- sigma0.seq[i]
  lambda0 <- lambda0.seq[j]
  
  assign("mu",mu,envir = .GlobalEnv)
  assign("sigma0.placebo",sigma0,envir = .GlobalEnv)
  assign("lambda0.placebo",lambda0,envir = .GlobalEnv)
  
  prior.placebo <- function(x){exp(-(abs(x-mu)/sigma0)^lambda0)}
  nc <- integrate(prior.placebo,lower=0+epsilon,upper=1-epsilon)[[1]]
  prior.placebo.nc <- function(x){exp(-(abs(x-mu)/sigma0)^lambda0)/nc}
  
  return(prior.placebo.nc)
}

skpt_tail_area <- function(){
  
  sigma0.seq <- seq(.01,0.5,by=0.01)
  lambda0.seq <- seq(1,2,by=0.1)
  result1 <- matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  result2 <- matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  q.outer <- 0.025  # y > x + delta.enth
  q.inner <- 0.010  # y > x + delta.intr
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.skpt <- sigma0.seq[i]
      lambda0.skpt <- lambda0.seq[j]
      
      prior.skpt <- function(x,y){
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
          exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)}
      nc <- tryCatch(integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
 error=function(e) integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
      
      result1[i,j] <- tryCatch(integral2(prior.skpt,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,singular=T)[[1]]/nc,
                             error=function(e) integral2(prior.skpt,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,abstol=1E-6)[[1]]/nc)
      result2[i,j] <- tryCatch(integral2(prior.skpt,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/nc, 
                             error=function(e) integral2(prior.skpt,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/nc)
    }
  }
  
  result3=abs(result1-q.outer)+abs(result2-q.inner)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  sigma0.skpt <- sigma0.seq[i]
  lambda0.skpt <- lambda0.seq[j]
  print(result1[i,j])
  print(result2[i,j])
  
  assign("sigma0.skpt",sigma0.skpt,envir = .GlobalEnv)
  assign("lambda0.skpt",lambda0.skpt,envir = .GlobalEnv)
  
  prior.skpt <- function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)}
  nc <- tryCatch(integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
               error=function(e) integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  prior.skpt.nc <- function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/nc}
  
}

enth_tail_area <- function(){
  
  sigma0.seq <- seq(0.085,0.90,by=0.001)
  lambda0.seq <- seq(2,2,by=0.1)
  result1 <- matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  #result2 <- matrix(NA,nrow=length(sigma0.seq),ncol=length(lambda0.seq))
  
  q.outer <- 0.025  # y > x + delta.enth
  #q.inner <- 0.010  # y > x + delta.intr
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.enth <- sigma0.seq[i]
      lambda0.enth <- lambda0.seq[j]
      
      prior.enth <- function(x,y){
        exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
          exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)}
      nc <- tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
 error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
      
      result1[i,j] <- tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]]/nc,
                             error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,abstol=1E-6)[[1]]/nc)      
      #result2[i,j] <- tryCatch(integral2(prior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/prior.enth.nc,
      #     error=function(e) integral2(prior.enth,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/prior.enth.nc)
    }
  }
  # 1-result1 since looking for q.outer area that y < x
  # 1-result2 since looking for q.inner area that y < x + delta.intr
  result3=abs(1-result1-q.outer)#+abs(1-result2-q.inner)
  index <- which(result3  ==  min(result3), arr.ind = TRUE)
  i <- index[1]
  j <- index[2]
  sigma0.enth <- sigma0.seq[i]
  lambda0.enth <- lambda0.seq[j]
  # 1-result1 since looking for q.outer area that y < x
  # 1-result2 since looking for q.inner area that y < x + delta.intr
  print(1-result1[i,j]-q.outer)
  #print(1-result2[i,j])
  
  assign("sigma0.enth",sigma0.enth,envir = .GlobalEnv)
  assign("lambda0.enth",lambda0.enth,envir = .GlobalEnv)
  
  prior.enth <- function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)}
  nc <- tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
  prior.enth.nc <- function(x,y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/nc}
  
  
  return(prior.enth.nc)
}