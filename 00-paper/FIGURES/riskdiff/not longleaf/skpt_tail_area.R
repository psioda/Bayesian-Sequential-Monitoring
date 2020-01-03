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
prior.skpt.nc<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
      error=function(e) integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])

result1[i,j]<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,singular=T)[[1]]/prior.skpt.nc,
     error=function(e) integral2(prior.skpt,xmin=0,xmax=1-delta.enth,ymin=function(x) x+delta.enth,ymax=1,abstol=1E-6)[[1]]/prior.skpt.nc)
result2[i,j]<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,singular=T)[[1]]/prior.skpt.nc, 
     error=function(e) integral2(prior.skpt,xmin=0,xmax=1-delta.intr,ymin=function(x) x+delta.intr,ymax=1,abstol=1E-6)[[1]]/prior.skpt.nc)
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
}