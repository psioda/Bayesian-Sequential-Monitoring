rm(list = ls())

## posterior coverage probability
## assuming equal mixture of priors
posterior.cov<-function(a1,b1,a2,b2,y0,y1,p,cred.tail){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
     beta(a2+y1,b2+y0)/beta(a2,b2))
  
  lower_cr<-0
  tail <-0
  while(tail<cred.tail/2){
    tail<-c*pbeta(lower_cr,a1+y1,b1+y0)+
      (1-c)*pbeta(lower_cr,a2+y1,b2+y0)
    if (tail<=cred.tail/2) lower_cr<-lower_cr+1e-3
  }
  
  upper_cr<-1
  tail <-0
  while(tail<cred.tail/2){
    tail<-c*pbeta(upper_cr,a1+y1,b1+y0,lower.tail=FALSE)+
      (1-c)*pbeta(upper_cr,a2+y1,b2+y0,lower.tail=FALSE)
    if (tail<=cred.tail/2) upper_cr<-upper_cr-1e-3
  }
  
  result<-(p>=lower_cr & p<=upper_cr)
  
  return(result)}
## posterior coverage probability
posterior.cov.2<-function(a.s.1,b.s.1,a.s.2,b.s.2,w.s.1,
                          a.e.1,b.e.1,a.e.2,b.e.2,w.e.1,
                          y0,y1,p,cred.tail){
  
  k.0<-1/2*c(w.s.1,1-w.s.1,w.e.1,1-w.e.1)
  
  c<-c(beta(a.s.1+y1,b.s.1+y0)/beta(a.s.1,b.s.1),
       beta(a.s.2+y1,b.s.2+y0)/beta(a.s.2,b.s.2),
       beta(a.e.1+y1,b.e.1+y0)/beta(a.e.1,b.e.1),
       beta(a.e.2+y1,b.e.2+y0)/beta(a.e.2,b.e.2))
  
  k.1<-k.0*c/(sum(k.0*c))
  
  lower_cr<-0
  tail <-0
  while(tail<cred.tail/2){
    tail<-k.1[1]*pbeta(lower_cr,a.s.1+y1,b.s.1+y0)+
          k.1[2]*pbeta(lower_cr,a.s.2+y1,b.s.2+y0)+
          k.1[3]*pbeta(lower_cr,a.e.1+y1,b.e.1+y0)+
          k.1[4]*pbeta(lower_cr,a.e.2+y1,b.e.2+y0)
    if (tail<=cred.tail/2) lower_cr<-lower_cr+1e-3
  }
  
  upper_cr<-1
  tail <-0
  while(tail<cred.tail/2){
    tail<-k.1[1]*pbeta(upper_cr,a.s.1+y1,b.s.1+y0,lower.tail=FALSE)+
          k.1[2]*pbeta(upper_cr,a.s.2+y1,b.s.2+y0,lower.tail=FALSE)+
          k.1[3]*pbeta(upper_cr,a.e.1+y1,b.e.1+y0,lower.tail=FALSE)+
          k.1[4]*pbeta(upper_cr,a.e.2+y1,b.e.2+y0,lower.tail=FALSE)
    if (tail<=cred.tail/2) upper_cr<-upper_cr-1e-3
  }
  
  result<-(p>=lower_cr & p<=upper_cr)

  return(result)
}

## posterior cdf
posterior.cdf<-function(a1,b1,a2,b2,y0,y1,q){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*pbeta(q,a1+y1,b1+y0)+(1-c)*pbeta(q,a2+y1,b2+y0)
  
  return(result)}

posterior.cdf.2<-function(a.s.1,b.s.1,a.s.2,b.s.2,w.s.1,
                          a.e.1,b.e.1,a.e.2,b.e.2,w.e.1,
                          y0,y1,q){
  
  k.0<-1/2*c(w.s.1,1-w.s.1,w.e.1,1-w.e.1)
  
  c<-c(beta(a.s.1+y1,b.s.1+y0)/beta(a.s.1,b.s.1),
       beta(a.s.2+y1,b.s.2+y0)/beta(a.s.2,b.s.2),
       beta(a.e.1+y1,b.e.1+y0)/beta(a.e.1,b.e.1),
       beta(a.e.2+y1,b.e.2+y0)/beta(a.e.2,b.e.2))
  
  k.1<-k.0*c/(sum(k.0*c))

  result<-k.1[1]*pbeta(q,a.s.1+y1,b.s.1+y0)+
          k.1[2]*pbeta(q,a.s.2+y1,b.s.2+y0)+
          k.1[3]*pbeta(q,a.e.1+y1,b.e.1+y0)+
          k.1[4]*pbeta(q,a.e.2+y1,b.e.2+y0)
  
  return(result)
}

## posterior mean
posterior.mean<-function(a1,b1,a2,b2,y0,y1){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*(a1+y1)/(a1+b1+y1+y0)+(1-c)*(a2+y1)/(a2+b2+y1+y0)
  
  return(result)}
posterior.mean.2<-function(a.s.1,b.s.1,a.s.2,b.s.2,w.s.1,
                           a.e.1,b.e.1,a.e.2,b.e.2,w.s.2,
                           y0,y1){
  
  k.0<-1/2*c(w.s.1,1-w.s.1,w.e.1,1-w.e.1)
  
  c<-c(beta(a.s.1+y1,b.s.1+y0)/beta(a.s.1,b.s.1),
       beta(a.s.2+y1,b.s.2+y0)/beta(a.s.2,b.s.2),
       beta(a.e.1+y1,b.e.1+y0)/beta(a.e.1,b.e.1),
       beta(a.e.2+y1,b.e.2+y0)/beta(a.e.2,b.e.2))
  
  k.1<-k.0*c/(sum(k.0*c))
  
  result<-k.1[1]*(a.s.1+y1)/(a.s.1+b.s.1+y1+y0)+
          k.1[2]*(a.s.2+y1)/(a.s.2+b.s.2+y1+y0)+
          k.1[3]*(a.e.1+y1)/(a.e.1+b.e.1+y1+y0)+
          k.1[4]*(a.e.2+y1)/(a.e.2+b.e.2+y1+y0)
  
  return(result)}