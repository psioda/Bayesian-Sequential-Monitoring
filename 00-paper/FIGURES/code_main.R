## Design parameters ##
p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
p.enth<-0.40
p.intr<-0.30
tail.skpt<-0.045  # tail probabilities for priors (low, high)
tail.enth<-0.05
sig.fut<-0.85     # significant trial result threshold
sig.eff<-0.95
cred.tail<-0.05   # credible interval is 1-cred.tail
max.ss<-76        # maximum sample size

## Simulation parameters ##
reps<-1000       # number of simulated trials per design
p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05) # range of response proportion
#p.range<-0.2                             
freq.mntr<-2      # frequency of monitoring
#freq.mntr<-c(1,2,4,8,16,76)
#freq.mntr<-rep(c(1,2,4,8,16,80),4) 
enr.shape<-1      # shape gamma dist enrollment
#enr.shape<-c(rep(1,12),rep(0.25,12))
#enr.shape<-rep(1,6)
out.mean<-4       # mean normal dist outcome
#out.mean<-rep(4,6)
#out.mean<-rep(c(rep(4,6),rep(8,6)),2)

spike<-0 # spike/slab version or regular version

#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,100,by=0.001)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.skpt
quantiles.skpt<-qbeta(tail.skpt,(p.skpt)*phi.range,(1-(p.skpt))*phi.range,
                      lower.tail=FALSE)
# lower tail probability equal to tail.enth
quantiles.enth<-qbeta(tail.enth,(p.enth)*phi.range,(1-(p.enth))*phi.range,
                      lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi.range[which.min(abs(p.enth-quantiles.skpt))] # fixed 5/13/19
phi_H<-phi.range[which.min(abs(p.skpt-quantiles.enth))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.skpt<-(p.skpt)*phi_L
beta.skpt<-(1-(p.skpt))*phi_L
alpha.enth<-(p.enth)*phi_H
beta.enth<-(1-(p.enth))*phi_H

#################################################################################################
## SIMULATIONS ##################################################################################
#################################################################################################

matrix.names<-c("eff.mon.initial","eff.mon.final",
                "eff.inf.initial","eff.inf.final",
                "fut.mon.initial","fut.mon.final",
                "fut.inf.initial","fut.inf.final",
                "phat.initial","phat.final",
                "ss.initial","ss.final",
                "post.mean.initial","post.mean.final",
                "cov.initial","cov.final")

outer<- array(NA,dim = c(length(freq.mntr),length(p.range),length(matrix.names)),
              dimnames = list(seq_len(length(freq.mntr)),p.range,matrix.names))

inner <- array(NA, dim=c(reps,length(matrix.names)), dimnames = list(seq_len(reps),matrix.names))

## posterior probability
probs.p<-c(0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1)
outer.p <- array(NA,dim=c(length(freq.mntr),length(p.range),length(probs.p)),
                 dimnames= list(seq_len(length(freq.mntr)),p.range,probs.p))
inner.p <- array(NA, dim=c(reps,2),dimnames=list(seq_len(reps),c("initial.p","final.p")))

outer.p.agree<-array(NA,dim=c(length(freq.mntr),length(p.range),3),
                     dimnames= list(seq_len(length(freq.mntr)),p.range,c("p.agree","efficacy","conditional")))

for (i in 1:length(freq.mntr)){

for (j in 1:length(p.range)){
  
  for (k in 1:reps){
    
    if (k%%1000==0){print(paste0("Model ",i,", Response p=",p.range[j],", Simulation ",k))}
    
    enr.times<-cumsum(rgamma(n=max.ss,shape=enr.shape[i],scale=0.5))
    outcome.times<-sort(enr.times+rnorm(n=max.ss,mean=out.mean[i],sd=0.25))
    
    responses<-rbinom(n=max.ss,size=1,prob=p.range[j])
    y1<-cumsum(responses)
    y0<-seq(1:length(responses))-y1
    
    ## monitoring
    if (spike==0){
    futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
    efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)}
    if (spike==1){
    futility<-posterior.quantile(a.s=alpha.enth.1,b.s=beta.enth.1,
                                 a.e=alpha.enth.2,b.e=beta.enth.2,
                                 y0=y0,y1=y1,q=p.intr,sims=sims)
      
    efficacy<-posterior.quantile(a.s=alpha.skpt.1,b.s=beta.skpt.1,
                                 a.e=alpha.skpt.2,b.e=beta.skpt.2,
                                 y0=y0,y1=y1,q=p.skpt,sims=sims)}

    n.initial<-min(
      which((futility>sig.fut) | ((1-efficacy)>sig.eff) & (seq(1:max.ss)%%freq.mntr==0)),
      max.ss,na.rm=TRUE)
    
    cutoff.time<-outcome.times[n.initial]
    n.final<-length(responses[enr.times<=cutoff.time])

    time<-c("initial","final")
    n<-c(n.initial,n.final)
    for (l in 1:2){
    inner[k,paste0("fut.mon.",time[l])]<-(futility[n[l]]>sig.fut)
    inner[k,paste0("eff.mon.",time[l])]<-(1-efficacy[n[l]]>sig.eff)
    inner[k,paste0("ss.",time[l])]<-n[l]
    inner[k,paste0("phat.",time[l])]<-y1[n[l]]/n[l]
    if (spike==0){
    inner[k,paste0("fut.inf.",time[l])]<-(posterior.quantile(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y0=y0[n[l]],y1=y1[n[l]],q=p.intr)>sig.fut)
    inner[k,paste0("eff.inf.",time[l])]<-(1-posterior.quantile(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y0=y0[n[l]],y1=y1[n[l]],q=p.skpt)>sig.eff)
    inner[k,paste0("post.mean.",time[l])]<-posterior.mean(a1=alpha.skpt,b1=beta.skpt,
                                                     a2=alpha.enth,b2=beta.enth,
                                                     y1=y1[n[l]],y0=y0[n[l]])
    inner[k,paste0("cov.",time[l])]<-posterior.cov(a1=alpha.skpt,b1=beta.skpt,
                                              a2=alpha.enth,b2=beta.enth,
                                              y1=y1[n[l]],y0=y0[n[l]],
                                              p=p.range[j],
                                              cred.tail=cred.tail)}
    if (spike==1){
      inner[k,paste0("post.mean.",time[l])]<-posterior.mean.2(
        a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
        a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,const.1=mix.1,
        a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
        a.e.2=alpha.enth.2,b.e.2=beta.enth.2,const.2=mix.2,
        y1=y1[n[l]],y0=y0[n[l]])
      inner.cov.initial[i]<-posterior.cov.2(
        a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
        a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,const.1=mix.1,
        a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
        a.e.2=alpha.enth.2,b.e.2=beta.enth.2,const.2=mix.2,
        y1=y1[n[l]],y0=y0[n[l]],p=p.range[j],cred.tail=cred.tail,sims=sims)}
    # posterior probability
    inner.p[k,paste0(time[l],".p")]<-(1-efficacy[n[l]])
    }
    
  }
outer[i,j,]<-apply(inner,MARGIN=2,FUN=mean)
outer.p[i,j,]<-quantile(inner.p[,"final.p"][inner.p[,"initial.p"]>sig.eff],probs=probs.p)
outer.p.agree[i,j,"p.agree"]<-sum((inner.p[,"initial.p"]>sig.eff)==(inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree[i,j,"efficacy"]<-sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree[i,j,"conditional"]<-sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/
                                                                            sum(inner.p[,"initial.p"]>sig.eff)
}
}