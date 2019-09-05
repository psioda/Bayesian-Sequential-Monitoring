

makeHist<-function(x){
  x<-seq(0,1,length=100)
  y<-rep(NA,length=length(x))
  for (i in 1:length(x)){
    if (x[i]< 0.5) y[i]=0.2
    if (x[i]>=0.5) y[i]=0.4
  }
  return(y)
}
  
makeHist(seq(0,1,by=0.01))
  
y
