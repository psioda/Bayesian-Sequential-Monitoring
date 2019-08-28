rm(list = ls())
## posterior coverage probability
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
posterior.cov.2<-function(a.s.1,b.s.1,a.s.2,b.s.2,const.1,
                          a.e.1,b.e.1,a.e.2,b.e.2,const.2,
                          y0,y1,p,cred.tail,sims){
  
  tmp1<-const.1*(beta(a.s.1+y1,b.s.1+y0)/beta(a.s.1,b.s.1))+
    (1-const.1)*(beta(a.s.2+y1,b.s.2+y0)/beta(a.s.2,b.s.2))
  tmp2<-const.2*(beta(a.e.1+y1,b.e.1+y0)/beta(a.e.1,b.e.1))+
    (1-const.2)*(beta(a.e.2+y1,b.e.2+y0)/beta(a.e.2,b.e.2))
  
  c<-(tmp1)/(tmp1+tmp2)
  
  theta.1.1<-rbeta(round(sims*c*const.1,log10(sims)),a.s.1+y1,b.s.1+y0)
  theta.1.2<-rbeta(round(sims*c*(1-const.1),log10(sims)),a.s.2+y1,b.s.2+y0)
  theta.2.1<-rbeta(round(sims*(1-c)*const.2,log10(sims)),a.e.1+y1,b.e.1+y0)
  theta.2.2<-rbeta(round(sims*(1-c)*(1-const.2),log10(sims)),a.e.2+y1,b.e.2+y0)
  theta<-sort(c(theta.1.1,theta.1.2,theta.2.1,theta.2.2))
  
  q.lower<-theta[round((cred.tail/2)*sims)]
  q.upper<-theta[round((1-cred.tail/2)*sims)]
  
  result<-(p>q.lower & p<q.upper)
  
  return(result)
}

## posterior quantile
posterior.quantile<-function(a1,b1,a2,b2,y0,y1,q){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
       beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*pbeta(q,a1+y1,b1+y0)+
    (1-c)*pbeta(q,a2+y1,b2+y0)
  
  return(result)}
posterior.quantile.2<-function(a.s.1,b.s.1,a.s.2,b.s.2,const.1,
                               a.e.1,b.e.1,a.e.2,b.e.2,const.2,
                               y0,y1,q,sims){
  
  tmp1<-const.1*(beta(a.s.1+y1,b.s.1+y0)/beta(a.s.1,b.s.1))+
    (1-const.1)*(beta(a.s.2+y1,b.s.2+y0)/beta(a.s.2,b.s.2))
  tmp2<-const.2*(beta(a.e.1+y1,b.e.1+y0)/beta(a.e.1,b.e.1))+
    (1-const.2)*(beta(a.e.2+y1,b.e.2+y0)/beta(a.e.2,b.e.2))
  
  c<-(tmp1)/(tmp1+tmp2)
  
  theta.1.1<-rbeta(round(sims*c*const.1,log10(sims)),a.s.1+y1,b.s.1+y0)
  theta.1.2<-rbeta(round(sims*c*(1-const.1),log10(sims)),a.s.2+y1,b.s.2+y0)
  theta.2.1<-rbeta(round(sims*(1-c)*const.2,log10(sims)),a.e.1+y1,b.e.1+y0)
  theta.2.2<-rbeta(round(sims*(1-c)*(1-const.2),log10(sims)),a.e.2+y1,b.e.2+y0)
  theta<-sort(c(theta.1.1,theta.1.2,theta.2.1,theta.2.2))
  
  result<-sum(theta<q)/sims
  
  return(result)
}

## posterior mean
posterior.mean<-function(a1,b1,a2,b2,y0,y1){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
       beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*(a1+y1)/(a1+b1+y1+y0)+
    (1-c)*(a2+y1)/(a2+b2+y1+y0)
  
  return(result)}
posterior.mean.2<-function(a.s.1,b.s.1,a.s.2,b.s.2,const.1,
                           a.e.1,b.e.1,a.e.2,b.e.2,const.2,
                           y0,y1){
  
  tmp1<-const.1*(beta(a.s.1+y1,b.s.1+y0)/beta(a.s.1,b.s.1))+
    (1-const.1)*(beta(a.s.2+y1,b.s.2+y0)/beta(a.s.2,b.s.2))
  tmp2<-const.2*(beta(a.e.1+y1,b.e.1+y0)/beta(a.e.1,b.e.1))+
    (1-const.2)*(beta(a.e.2+y1,b.e.2+y0)/beta(a.e.2,b.e.2))
  
  c<-(tmp1)/(tmp1+tmp2)
  
  result<-c*(const.1*(a.s.1+y1)/(a.s.1+b.s.1+y1+y0)+
               (1-const.1)*(a.s.2+y1)/(a.s.2+b.s.2+y1+y0))+
    (1-c)*(const.2*(a.e.1+y1)/(a.e.1+b.e.1+y1+y0)+
             (1-const.2)*(a.e.2+y1)/(a.e.2+b.e.2+y1+y0))
  
  return(result)}