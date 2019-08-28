## posterior mean
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


## posterior coverage probability
## inference prior is a mixture of skeptical and enthuastic,
## which themselves are mixtures of beta distributions
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

# alpha.skpt.1<-alpha.skpt
# beta.skpt.1<-beta.skpt
# alpha.skpt.2<-alpha.skpt
# beta.skpt.2<-beta.skpt
# mix.1<-.5
# alpha.enth.1<-alpha.enth
# beta.enth.1<-beta.enth
# alpha.enth.2<-alpha.enth
# beta.enth.2<-beta.enth
# mix.2<-.5
# 
# alpha.skpt.1<-1.59
# beta.skpt.1<-6.36
# alpha.skpt.2<-20
# beta.skpt.2<-80
# mix.1<-0.5
# alpha.enth.1<-3.56
# beta.enth.1<-5.34
# alpha.enth.2<-12
# beta.enth.2<-18
# mix.2<-0.5

# alpha.skpt.1<-1.7752
# beta.skpt.1<-10.9048
# alpha.skpt.2<-5.7226
# beta.skpt.2<-16.2874
# mix.1<-0.5
# alpha.enth.1<-8.3979
# beta.enth.1<-18.6921
# alpha.enth.2<-11.4219
# beta.enth.2<-11.8881
# mix.2<-0.5

# futility<-posterior.quantile(a.s=alpha.enth.1,b.s=beta.enth.1,
#                              a.e=alpha.enth.2,b.e=beta.enth.2,
#                              y0=y0,y1=y1,q=p.intr,sims=sims)

#efficacy<-posterior.quantile(a.s=alpha.skpt.1,b.s=beta.skpt.1,
#                             a.e=alpha.skpt.2,b.e=beta.skpt.2,
#                             y0=y0,y1=y1,q=p.skpt,sims=sims)

# 8/2/19

# inner.fut.final[i]<-(posterior.quantile(a.s=alpha.enth.1,b.s=beta.enth.1,
#                                         a.e=alpha.enth.2,b.e=beta.enth.2,
#                                         y0=y0.final,y1=y1.final,
#                                         q=p.intr,sims=sims)
#                      >sig.fut)
# inner.eff.final[i]<-((1-posterior.quantile(a.s=alpha.skpt.1,b.s=beta.skpt.1,
#                                            a.e=alpha.skpt.2,b.e=beta.skpt.2,
#                                            y0=y0.final,y1=y1.final,
#                                            q=p.skpt,sims=sims))
#                      >sig.eff)





# inner.eff.final.mix[i]<-((1-posterior.quantile(a.s=alpha.skpt,b.s=beta.skpt,
#                                            a.e=alpha.enth,b.e=beta.enth,
#                                            y0=y0.final,y1=y1.final,
#                                            q=p.skpt,sims=sims))
#                      >sig.eff)


# posterior mean
# inner.post.mean.initial[i]<-posterior.mean.2(
#   a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
#   a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,const.1=mix.1,
#   a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
#   a.e.2=alpha.enth.2,b.e.2=beta.enth.2,const.2=mix.2,
#   y1=y1,y0=y0)
# 
# inner.post.mean.final[i]<-posterior.mean.2(
#   a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
#   a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,const.1=mix.1,
#   a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
#   a.e.2=alpha.enth.2,b.e.2=beta.enth.2,const.2=mix.2,
#   y1=y1.final,y0=y0.final)
# 
# inner.cov.initial[i]<-posterior.cov.2(
#   a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
#   a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,const.1=mix.1,
#   a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
#   a.e.2=alpha.enth.2,b.e.2=beta.enth.2,const.2=mix.2,
#   y1=y1,y0=y0,p=p.range[j],cred.tail=cred.tail,sims=sims)
# 
# inner.cov.final[i]<-posterior.cov.2(
#   a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
#   a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,const.1=mix.1,
#   a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
#   a.e.2=alpha.enth.2,b.e.2=beta.enth.2,const.2=mix.2,
#   y1=y1.final,y0=y0.final,p=p.range[j],cred.tail=cred.tail,sims=sims)


## posterior quantile
posterior.quantile<-function(a1,b1,a2,b2,y0,y1,q){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
       beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*pbeta(q,a1+y1,b1+y0)+
    (1-c)*pbeta(q,a2+y1,b2+y0)
  
  return(result)}


futility.inf<-poste


responses<-c(rep(c(1,0,1,0),5))
y1<-cumsum(responses)
y0<-seq(1:length(responses))-y1

futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
futility.inf<-posterior.quantile(a1=alpha.skpt,b1=beta.skpt,
                                 a2=alpha.enth,b2=beta.enth,
                                 y0=y0,y1=y1,q=p.intr)
efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)
efficacy.inf<-posterior.quantile(a1=alpha.skpt,b1=beta.skpt,
                                 a2=alpha.enth,b2=beta.enth,
                                 y0=y0,y1=y1,q=p.skpt)




if ((n%%freq.mntr[k]==0 | n==max.ss) & futility>sig.fut){
  inner[i,"fut"]<-1
  inner[i,"eff"]<-0
  inner[i,"inc"]<-0
  break
}
else if ((n%%freq.mntr[k]==0 | n==max.ss) & (1-efficacy)>sig.eff){
  inner[i,"eff"]<-1
  inner[i,"fut"]<-0
  inner[i,"inc"]<-0
  break
}
else {
  inner[i,"inc"]<-1
  inner[i,"eff"]<-0
  inner[i,"fut"]<-0
}
}

cutoff.time<-final.times[n]
#responses.final<-responses[event.times<=cutoff.time]
responses.final<-responses[1:n]
n.final<-length(responses.final)
y1.final<-sum(responses.final)
y0.final<-n.final-y1.final

inner[i,"fut.final"]<-(pbeta(p.intr,alpha.enth+y1.final,beta.enth+y0.final,
                             lower.tail=TRUE)>sig.fut)
inner[i,"eff.final"]<-(pbeta(p.skpt,alpha.skpt+y1.final,beta.skpt+y0.final,
                             lower.tail=FALSE)>sig.eff)
inner[i,"inc.final"]<-1-inner[i,"fut.final"]-inner[i,"eff.final"]

inner[i,"ss.initial"]<-n
inner[i,"phat.initial"]<-y1/n
inner[i,"ss.final"]<-n.final
inner[i,"phat.final"]<-y1.final/n.final


inner[i,"post.mean.initial"]<-posterior.mean(a1=alpha.skpt,b1=beta.skpt,
                                             a2=alpha.enth,b2=beta.enth,
                                             y1=y1,y0=y0)

inner[i,"post.mean.final"]<-posterior.mean(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y1=y1.final,y0=y0.final)

inner[i,"cov.initial"]<-posterior.cov(a1=alpha.skpt,b1=beta.skpt,
                                      a2=alpha.enth,b2=beta.enth,
                                      y1=y1,y0=y0,
                                      p=p.range[j],
                                      cred.tail=cred.tail)

inner[i,"cov.final"]<-posterior.cov(a1=alpha.skpt,b1=beta.skpt,
                                    a2=alpha.enth,b2=beta.enth,
                                    y1=y1.final,y0=y0.final,
                                    p=p.range[j],
                                    cred.tail=cred.tail)



