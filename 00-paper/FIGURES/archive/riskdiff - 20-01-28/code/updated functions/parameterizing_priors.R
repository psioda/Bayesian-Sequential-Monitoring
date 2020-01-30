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