rm(list = ls())

# mean of skeptical prior
p.skpt<-0.20
# mean of enthuastic prior
p.enth<-0.40
# futility theta
p.intr<-0.30
# value of true response proportion
p.range<-0.2 # only interested in type 1 error
#p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05)
#p.range<-seq(p.skpt,p.enth,by=0.1)
# tail probabilities for priors (low, high)
tail.skpt<-0.045
tail.enth<-0.05
# maximum sample sizes
max.ss<-80
# significant trial result threshold
sig.fut<-0.85
sig.eff<-0.95
# number of simulated trials per design
reps<-50000
# compute empirical quantities for credible interval
sims<-10000
# credible interval is 1-cred.tail
cred.tail<-0.05
# design parameters
# frequency of monitoring
#freq.mntr<-rep(2)
freq.mntr<-rep(c(1,2,4,8,16,80),4) 
# shape gamma dist enrollment
#enr.shape<-1
enr.shape<-c(rep(1,12),rep(0.25,12))
# mean normal dist outcome
#out.mean<-4
out.mean<-rep(c(rep(4,6),rep(8,6)),2)


## posterior coverage probability
## inference prior is a mixture of skeptical and enthuastic
posterior.cov<-function(a.s,b.s,a.e,b.e,y0,y1,p,cred.tail,sims){
  
  c<-beta(a.s+y1,b.s+y0)/beta(a.s,b.s)/
    (beta(a.s+y1,b.s+y0)/beta(a.s,b.s)+
     beta(a.e+y1,b.e+y0)/beta(a.e,b.e))
  
  theta.1<-rbeta(round(sims*c,log10(sims)),a.s+y1,b.s+y0)
  theta.2<-rbeta(round(sims*(1-c),log10(sims)),a.e+y1,b.e+y0)
  theta<-sort(c(theta.1,theta.2))
  
  q.lower<-theta[round((cred.tail/2)*sims)]
  q.upper<-theta[round((1-cred.tail/2)*sims)]

  result<-(p>q.lower & p<q.upper)

  return(result)}

## posterior quantile
## inference prior is a mixture of skeptical and enthuastic
posterior.quantile<-function(a.s,b.s,a.e,b.e,y0,y1,q,sims){
  
  c<-beta(a.s+y1,b.s+y0)/beta(a.s,b.s)/
    (beta(a.s+y1,b.s+y0)/beta(a.s,b.s)+
     beta(a.e+y1,b.e+y0)/beta(a.e,b.e))
  
  theta.1<-rbeta(round(sims*c,log10(sims)),a.s+y1,b.s+y0)
  theta.2<-rbeta(round(sims*(1-c),log10(sims)),a.e+y1,b.e+y0)
  theta<-sort(c(theta.1,theta.2))
  
  result<-sum(theta<q)/sims
  
  return(result)}

## posterior mean
## inference prior is a mixture of skeptical and enthuastic
posterior.mean<-function(a.s,b.s,a.e,b.e,y0,y1){
  
  c<-beta(a.s+y1,b.s+y0)/beta(a.s,b.s)/
     (beta(a.s+y1,b.s+y0)/beta(a.s,b.s)+
      beta(a.e+y1,b.e+y0)/beta(a.e,b.e))
  
  result<-c*(a.s+y1)/(a.s+b.s+y1+y0)+
      (1-c)*(a.e+y1)/(a.e+b.e+y1+y0)
  
  return(result)}

## posterior mean
## inference prior is a mixture of skeptical and enthuastic,
## which themselves are mixtures of beta distributions
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
## inference prior is a mixture of skeptical and enthuastic,
## which themselves are mixtures of beta distributions
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
#example1<-function(){

#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,100,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.skpt
quantiles.skpt<-qbeta(tail.skpt,(p.skpt)*phi.range,(1-(p.skpt))*phi.range,lower.tail=FALSE)
# lower tail probability equal to tail.enth
quantiles.enth<-qbeta(tail.enth,(p.enth)*phi.range,(1-(p.enth))*phi.range,lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi.range[which.min(abs(p.enth-quantiles.skpt))] # fixed 5/13/19
phi_H<-phi.range[which.min(abs(p.skpt-quantiles.enth))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.skpt<-(p.skpt)*phi_L
beta.skpt<-(1-(p.skpt))*phi_L
alpha.enth<-(p.enth)*phi_H
beta.enth<-(1-(p.enth))*phi_H

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

#################################################################################################
## SIMULATIONS ##################################################################################
#################################################################################################

## Step 1: Create outer loop based on frequency of interim analyses
# stop early for efficacy
outer.eff<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.eff.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.eff.final.mix<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
# stop early for futility
outer.fut<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.fut.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
# inconclusive findings
outer.inc<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.inc.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
# sample mean
outer.phat.initial<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.phat.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
# sample size
outer.ss.initial<-matrix(nrow=length(freq.mntr),ncol=length(p.range)) 
outer.ss.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range)) 
# posterior mean
outer.post.mean.initial<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.post.mean.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
# coverage probability
outer.cov.initial<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.cov.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))

for (k in 1:length(freq.mntr)){ # frequency of monitoring

for (j in 1:length(p.range)){ # true response proportion
  
  inner.eff<-vector(length=reps)
  inner.eff.final<-vector(length=reps)
  inner.eff.final.mix<-vector(length=reps)
  inner.fut<-vector(length=reps)
  inner.fut.final<-vector(length=reps)
  inner.inc<-vector(length=reps)
  inner.inc.final<-vector(length=reps)
  inner.ss.initial<-vector(length=reps)
  inner.ss.final<-vector(length=reps)
  inner.phat.initial<-vector(length=reps)
  inner.phat.final<-vector(length=reps)
  inner.post.mean.initial<-vector(length=reps)
  inner.post.mean.final<-vector(length=reps)
  inner.cov.initial<-vector(length=reps)
  inner.cov.final<-vector(length=reps)

  for (i in 1:reps){ # simulate specified trial design
    
    if (i%%1000==0){
      print(paste0("Model ",k)) 
      print(paste0("Response p=",p.range[j]))
      print(paste0("Simulation ",i))
      }
    
    efficacy<-0 # necessary?
    futility<-0 # necessary?
    inner.inc[i]<-1
    cutoff.time<-vector()
    event.times<-cumsum(rgamma(n=max.ss,shape=enr.shape[k],scale=0.5))
    final.times<-event.times+rnorm(n=max.ss,mean=out.mean[k],sd=0.25)
    responses<-rbinom(n=max.ss,size=1,prob=p.range[j])
    
    for (n in 1:max.ss){
      
      y1<-sum(responses[1:n])
      y0=n-y1
      
      # futility (even optimist would give up)
      futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
      # futility<-posterior.quantile(a.s=alpha.enth.1,b.s=beta.enth.1,
      #                              a.e=alpha.enth.2,b.e=beta.enth.2,
      #                              y0=y0,y1=y1,q=p.intr,sims=sims)
      
      # efficacy (even pessimist would accept)
      efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)
      #efficacy<-posterior.quantile(a.s=alpha.skpt.1,b.s=beta.skpt.1,
      #                             a.e=alpha.skpt.2,b.e=beta.skpt.2,
      #                             y0=y0,y1=y1,q=p.skpt,sims=sims)

      if ((n%%freq.mntr[k]==0 | n==max.ss) & futility>sig.fut){
        inner.fut[i]<-1
        inner.inc[i]<-0
        break
        }
      else if ((n%%freq.mntr[k]==0 | n==max.ss) & (1-efficacy)>sig.eff){
        inner.eff[i]<-1
        inner.inc[i]<-0
        break
        }
      else {
        inner.inc[i]<-1
      }
    }
    cutoff.time<-final.times[n]
    responses.final<-responses[event.times<=cutoff.time]
    n.final<-length(responses.final)
    y1.final<-sum(responses.final)
    y0.final<-n.final-y1.final
    
    inner.fut.final[i]<-(pbeta(p.intr,alpha.enth+y1.final,beta.enth+y0.final,
                              lower.tail=TRUE)>sig.fut)
    inner.eff.final[i]<-(pbeta(p.skpt,alpha.skpt+y1.final,beta.skpt+y0.final,
                              lower.tail=FALSE)>sig.eff)
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
    inner.inc.final[i]<-1-inner.fut.final[i]-inner.eff.final[i]
    
    
    
    
    # inner.eff.final.mix[i]<-((1-posterior.quantile(a.s=alpha.skpt,b.s=beta.skpt,
    #                                            a.e=alpha.enth,b.e=beta.enth,
    #                                            y0=y0.final,y1=y1.final,
    #                                            q=p.skpt,sims=sims))
    #                      >sig.eff)

    
    
    inner.ss.initial[i]<-n
    inner.phat.initial[i]<-y1/n
    inner.ss.final[i]<-n.final
    inner.phat.final[i]<-y1.final/n.final
    
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
    
    inner.post.mean.initial[i]<-posterior.mean(a.s=alpha.skpt,b.s=beta.skpt,
                                               a.e=alpha.enth,b.e=beta.enth,
                                               y1=y1,y0=y0)

    inner.post.mean.final[i]<-posterior.mean(a.s=alpha.skpt,b.s=beta.skpt,
                                             a.e=alpha.enth,b.e=beta.enth,
                                             y1=y1.final,y0=y0.final)

    # inner.cov.initial[i]<-posterior.cov(a.s=alpha.skpt,b.s=beta.skpt,
    #                                     a.e=alpha.enth,b.e=beta.enth,
    #                                     y1=y1,y0=y0,
    #                                     p=p.range[j],
    #                                     cred.tail=cred.tail,sims=sims)
    # 
    # inner.cov.final[i]<-posterior.cov(a.s=alpha.skpt,b.s=beta.skpt,
    #                                   a.e=alpha.enth,b.e=beta.enth,
    #                                   y1=y1.final,y0=y0.final,
    #                                   p=p.range[j],
    #                                   cred.tail=cred.tail,sims=sims)
  }

  outer.ss.initial[k,j]<-mean(inner.ss.initial)
  outer.ss.final[k,j]<-mean(inner.ss.final)
  outer.phat.initial[k,j]<-mean(inner.phat.initial)
  outer.phat.final[k,j]<-mean(inner.phat.final)
  outer.fut[k,j]<-mean(inner.fut)
  outer.fut.final[k,j]<-mean(inner.fut.final)
  outer.eff[k,j]<-mean(inner.eff)
  outer.eff.final[k,j]<-mean(inner.eff.final)
  outer.eff.final.mix[k,j]<-mean(inner.eff.final.mix)
  outer.inc[k,j]<-mean(inner.inc)
  outer.inc.final[k,j]<-mean(inner.inc.final)
  outer.post.mean.initial[k,j]<-mean(inner.post.mean.initial)
  outer.post.mean.final[k,j]<-mean(inner.post.mean.final)
  outer.cov.initial[k,j]<-mean(inner.cov.initial)
  outer.cov.final[k,j]<-mean(inner.cov.final)
}
}

par(mfrow = c(1,1)) 
k<-1
plot(p.range,outer.eff[k,],type='l',ylim=c(0,1),lwd=2,lty='longdash',
     ylab="Probability",xlab="",
     main="Probability of Efficacy\n 2 subjects per month, 4 month follow-up, analyze after every 2 results",
     axes=FALSE)
box()
text(p.range,outer.eff[k,],
     labels=format(round(outer.eff[k,],digits=2),nsmall=2),pos=3)
#lines(p.range,outer.fut.final[k,],col='red',lwd=2)
lines(p.range,outer.fut[k,],lwd=2,lty='longdash')
text(p.range,outer.fut[k,],
     labels=format(round(outer.fut[k,],digits=2),nsmall=2),pos=1)
#lines(p.range,outer.eff.final[k,],lwd=2,lty='dotted')
#lines(p.range,outer.eff.final.mix[k,],lwd=2)
lines(p.range,outer.inc[k,],lwd=2,lty='longdash')
text(p.range,outer.inc[k,],
     labels=format(round(outer.inc[k,],digits=2),nsmall=2),pos=1)
#lines(p.range,outer.inc.final[k,],lwd=2)
axis(1,las=0,at=p.range,labels=format(p.range,nsmall=2))
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')
abline(v=c(p.skpt,p.enth),col='grey',lty='dashed')
row<-1
  for (column in 1:length(p.range)){
    mtext(text=paste0(format(round(outer.ss.initial[k,column],digits=1),nsmall=1),
                      " + ",
                      format(round(outer.ss.final[k,column]-
                                  outer.ss.initial[k,column],digits=1),nsmall=1),
                      " = ",
                      format(round(outer.ss.final[k,column],digits=1),nsmall=1)),
                      side=1,line=row+1,at=p.range[column])
  }
row<-2
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",
        format(round(outer.post.mean.initial[k,column],digits=3),nsmall=3),
                    " (F) ",
        format(round(outer.post.mean.final[k,column],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[column])
}
row<-3
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",
        format(round(outer.cov.initial[k,column],digits=3),nsmall=3),
                    " (F) ",
        format(round(outer.cov.final[k,column],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[column])
}
column<-1
#legend(x=0.325,y=0.6,legend=c("Skeptic at Interim","Skeptical at Final","Mixture Prior"),
#       lty=c('longdash','dotted','solid'),
#       cex=0.8,
#       box.lwd = 1,box.col = "black",bg = "white",pt.cex = 1)
mtext(text="SS",side=1,line=2,at=0.125)
mtext(text="PM",side=1,line=3,at=0.125)
mtext(text="CP",side=1,line=4,at=0.125)
#### END EXAMPLES FOR 7/26/19 ####




#### Type I error graphs ####
par(mfrow = c(2, 2)) 
k1<-4
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),col='red',lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2])
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                   outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=8-column,
    at=1/freq.mntr[column])
}
k1<-8
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),col='red',lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2])
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                   outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=8-column,
    at=0.5)
}
k1<-4
k2<-0.25
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),col='red',lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2])
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
row<-1
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
                    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                 outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
                    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
      digits=1),nsmall=1)),
        side=1,
    line=8-column,
    at=0.5)
}
k1<-8
k2<-0.25
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),col='red',lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2])
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column]-
                   outer.ss.initial[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer.ss.final[out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=8-column,
    at=0.5)
}






outer_trial_result_binary # probability of rejecting null hypothesis (alpha at p.skpt, 1-beta at p.enth)
outer_trial_result
outer.ss.initial

# Step 5: Plot Skeptical and enthuastic priors separately
## 7/19/19 match slide 14 from FDA presentation, shade in gray scale, set 350 DPI
par(ask=TRUE)
par(mfrow = c(1, 2)) 
x<-seq(0,1,by=0.01)
# Make 2 boxplots
#low (skeptical)
plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",
     xlab="Response Probability",ylab="Density Value",
     main="Skeptical Beta Prior",
     #xaxt="n",
     ylim=c(0,max(dbeta(x,alpha.skpt,beta.skpt))))

polygon(c(x[x<=0.4],0.4),c(dbeta(x,alpha.skpt,beta.skpt)[x<=0.4],0),col="blue")
polygon(c(x[x>=0.4],0.4),c(dbeta(x,alpha.skpt,beta.skpt)[x>=0.4],0),col="red")
text(0.5, 4,expression(P(theta>0.40)==0.045),pos=4)
text(0.5,3.8,"Approximately N=14.0 Subjects",pos=4)
text(0.5,3.6,expression(E(theta)==0.20),pos=4)
#axis(1,at=p.skpt,labels=expression(theta[0]))
#axis(1,at=seq(0,1,by=0.2))
#high (enthuastic)
plot(x,dbeta(x,alpha.enth,beta.enth),type="l",
     xlab="Response Probability",
     ylab="",
     main="Enthuastic Beta Prior",
     #xaxt="n",
     yaxt="n",
     ylim=c(0,max(dbeta(x,alpha.enth,beta.enth))))
axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
polygon(c(x[x<=0.2],0.2),c(dbeta(x,alpha.enth,beta.enth)[x<=0.2],0),col="red")
polygon(c(x[x>=0.2],0.2),c(dbeta(x,alpha.enth,beta.enth)[x>=0.2],0),col="blue")
text(0.5, 3,expression(P(theta<0.20)==0.05),pos=4)
text(0.5,2.85,"Approximately N=14.0 Subjects",pos=4)
text(0.5,2.7,expression(E(theta)==0.40),pos=4)



# Step 6: Plot inference priors (mixtures)
## 7/19/19 make 50:50 solid black line, 25:75 dark grey with different dash patterns,
# make skeptical/enthuastic light grey with smaller thinkness
# no legend, put omega on actual lines

par(mfrow = c(1, 1)) 
plot(x,1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth),
     type="l",xlab="Response Probability",ylab="Density Value",
     main="50/50 Mixture of Skeptical and Enthuastic Priors",
     ylim=c(0,max(1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth))))
y<-1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth)
polygon(c(x,1),c(y,0),col="blue")
polygon(c(0.2,x[x>=0.2 & x<=0.4],0.4),c(0,y[x>=0.2 & x<=0.4],0),col="darkblue")


par(mfrow = c(1, 1)) 
plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",col="red",xlab="Response Probability",ylab="Density Value",main="Inference Priors",
     xaxt="n",
     ylim=c(0,max(dbeta(x,alpha.skpt,beta.skpt),dbeta(x,alpha.enth,beta.enth))))
axis(1,at=p.skpt,labels=expression(theta[0]))
lines(x,dbeta(x,alpha.enth,beta.enth),type="l",col="blue")
axis(1,at=p.enth,labels=expression(theta[A]))
axis(1,at=p.skpt,labels=expression(theta[0]))
points(x,1/4*dbeta(x,alpha.skpt,beta.skpt)+3/4*dbeta(x,alpha.enth,beta.enth),type='l',lty=3)
points(x,1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth),type='l',lty=2)
points(x,3/4*dbeta(x,alpha.skpt,beta.skpt)+1/4*dbeta(x,alpha.enth,beta.enth),type='l',lty=1)
legend(x=0.55,y=4.5,
       legend=c("Skeptical","75:25","50:50","25:75","Enthuastic"),
       col=c("red","black","black","black","blue"),lty=c(1,1,2,3,1),
       cex=0.8,text.font=3)

## Plot of power curves
par(ask=TRUE)
plot(p.range,outer_trial_result_binary[1,],type='l',
     xlab="Response Probability",ylab="Probability of Rejection",main="Power Curves")
lines(p.range,outer_trial_result_binary[2,],type='l')
lines(p.range,outer_trial_result_binary[3,],type='l')
lines(p.range,outer_trial_result_binary[4,],type='l')
legend(x=p.enth*(2/3),y=0.8,
       title="# enrolled",
       legend=freq.mntr,
       col=freq.mntr,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=TRUE)
plot(p.range,outer.ss.initial[4,],type='l',
     xlab="Response Probability",ylab="Expected Sample Size",main="Sample Size",
     ylim=c(max.ss/2,max.ss*1.1))
axis(2,at=max.ss,labels="Max")
points(p.range,outer.ss.initial[3,],type='l',lty=2)
points(p.range,outer.ss.initial[2,],type='l',lty=3)
points(p.range,outer.ss.initial[1,],type='l',lty=4)
legend(x=p.skpt*(3/2),max.ss*(4/5),
       title="# enrolled",
       legend=freq.mntr,
       col=freq.mntr,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=TRUE)
plot(p.range,abs(outer.phat.initial[4,]-p.range),type='l',
     xlab="Response Probability",ylab="Bias",main="Bias",
     ylim=c(0,0.1))
points(p.range,abs(outer.phat.initial[3,]-p.range),type='l',lty=2)
points(p.range,abs(outer.phat.initial[2,]-p.range),type='l',lty=3)
points(p.range,abs(outer.phat.initial[1,]-p.range),type='l',lty=4)
legend(x=p.enth*(6/7),0.1,
       title="# enrolled",
       legend=freq.mntr,
       col=freq.mntr,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=FALSE)

#}

#example1()
