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
      prior.enth.nc<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
            error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
      
      result1[i,j]<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,singular=T)[[1]]/prior.enth.nc,
           error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=function(x) x,ymax=1,abstol=1E-6)[[1]]/prior.enth.nc)      
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
}