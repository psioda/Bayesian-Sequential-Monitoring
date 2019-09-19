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

outer.p.agree<-array(NA,dim=c(length(freq.mntr),length(p.range),length(c("p.agree","efficacy","conditional"))),
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
    
    if (spike==0){
    futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
    efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)
    }
    
    if (spike==1){
    futility<-rep(NA,length=max.ss)
    efficacy<-rep(NA,length=max.ss)
    for (h in 1:max.ss){
      # skeptical posterior
      posterior.skpt<-function(x){
        exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
      }
      nc.skpt<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
      posterior.nc.skpt<-function(x){
        exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/nc.skpt
      }
      efficacy[h]<-integrate(posterior.nc.skpt,lower=0+epsilon,upper=p.skpt)[[1]]
      
      # enthuastic posterior
      posterior.enth<-function(x){
        exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
      }
      nc.enth<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
      posterior.nc.enth<-function(x){
        exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/nc.enth
      }
      futility[h]<-integrate(posterior.nc.enth,lower=0+epsilon,upper=p.intr)[[1]]
    }
    }
    
    n.initial<-min(
      which(((futility>sig.fut) | ((1-efficacy)>sig.eff)) & (seq(1:max.ss)%%freq.mntr[i]==0)),
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
    inner[k,paste0("fut.inf.",time[l])]<-(posterior.cdf(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y0=y0[n[l]],y1=y1[n[l]],q=p.intr,w=0.5)>sig.fut)
    inner[k,paste0("eff.inf.",time[l])]<-(1-posterior.cdf(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y0=y0[n[l]],y1=y1[n[l]],q=p.skpt,w=0.5)>sig.eff)
    inner[k,paste0("post.mean.",time[l])]<-posterior.mean(a1=alpha.skpt,b1=beta.skpt,
                                                     a2=alpha.enth,b2=beta.enth,
                                                     y1=y1[n[l]],y0=y0[n[l]])
    inner[k,paste0("cov.",time[l])]<-posterior.cov(a1=alpha.skpt,b1=beta.skpt,
                                              a2=alpha.enth,b2=beta.enth,
                                              y1=y1[n[l]],y0=y0[n[l]],
                                              p=p.range[j],
                                              cred.tail=cred.tail)
      }
    if (spike==1){
      ######################
      ### SKEPTICAL CASE ###
      ######################
      prior.skpt<-function(x){
        exp(-(abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt)
      }
      prior.skpt.nc<-integrate(prior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
      
      posterior.skpt<-function(x){
        exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
      }
      posterior.skpt.nc<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
      
      post.prob.skpt<-posterior.skpt.nc/prior.skpt.nc
      
      posterior.skpt.nc.exp<-function(x){
        x*exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/posterior.skpt.nc
      }
      skpt.exp<-integrate(posterior.skpt.nc.exp,lower=0+epsilon,upper=1-epsilon)[[1]]
      
      #######################
      ### ENTHUASTIC CASE ###
      #######################
      prior.enth<-function(x){
        exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
      }
      prior.enth.nc<-integrate(prior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
      
      posterior.enth<-function(x){
        exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
      }
      posterior.enth.nc<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
      
      post.prob.enth<-posterior.enth.nc/prior.enth.nc
      
      posterior.enth.nc.exp<-function(x){
        x*exp(y1[n[l]]*log(x)+y0[n[l]]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/posterior.enth.nc
      }
      enth.exp<-integrate(posterior.enth.nc.exp,lower=0+epsilon,upper=1-epsilon)[[1]]
      
      ####################
      ### WEIGHTED AVG ###
      ####################
      c<-post.prob.skpt/(post.prob.skpt+post.prob.enth)
      inner[k,paste0("post.mean.",time[l])]<-c*skpt.exp+(1-c)*enth.exp
    }
    # inner[k,paste0("fut.inf.",time[l])]<-(posterior.cdf.2(
    #                                       a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
    #                                       a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,
    #                                       w.s.1=mix.1,
    #                                       a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
    #                                       a.e.2=alpha.enth.2,b.e.2=beta.enth.2,
    #                                       w.e.1=mix.2,
    #                                       y1=y1[n[l]],y0=y0[n[l]],q=p.intr)>sig.fut)
    # 
    # inner[k,paste0("eff.inf.",time[l])]<-(1-posterior.cdf.2(
    #                                         a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
    #                                         a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,
    #                                         w.s.1=mix.1,
    #                                         a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
    #                                         a.e.2=alpha.enth.2,b.e.2=beta.enth.2,
    #                                         w.e.1=mix.2,
    #                                         y1=y1[n[l]],y0=y0[n[l]],q=p.skpt)>sig.eff)
    # 
    # inner[k,paste0("post.mean.",time[l])]<-posterior.mean.2(
    #     a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
    #     a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,w.s.1=mix.1,
    #     a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
    #     a.e.2=alpha.enth.2,b.e.2=beta.enth.2,w.e.1=mix.2,
    #     y1=y1[n[l]],y0=y0[n[l]])
    # 
    # inner[k,paste0("cov.",time[l])]<-posterior.cov.2(
    #     a.s.1=alpha.skpt.1,b.s.1=beta.skpt.1,
    #     a.s.2=alpha.skpt.2,b.s.2=beta.skpt.2,w.s.1=mix.1,
    #     a.e.1=alpha.enth.1,b.e.1=beta.enth.1,
    #     a.e.2=alpha.enth.2,b.e.2=beta.enth.2,w.e.1=mix.2,
    #     y1=y1[n[l]],y0=y0[n[l]],p=p.range[j],cred.tail=cred.tail)

    # posterior probability
    inner.p[k,paste0(time[l],".p")]<-(1-efficacy[n[l]])
    }
    
  }
outer[i,j,]<-apply(inner,MARGIN=2,FUN=mean)
outer.p[i,j,]<-quantile(inner.p[,"final.p"][inner.p[,"initial.p"]>sig.eff & inner.p[,"final.p"]<sig.eff],
                        probs=probs.p)
outer.p.agree[i,j,"p.agree"]<-sum((inner.p[,"initial.p"]>sig.eff)==(inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree[i,j,"efficacy"]<-sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree[i,j,"conditional"]<-sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/
                                                                            sum(inner.p[,"initial.p"]>sig.eff)
}
}