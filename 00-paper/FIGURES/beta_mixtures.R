rm(list = ls())

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
reps<-5000       # number of simulated trials per design
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

## posterior quantile
posterior.quantile<-function(a1,b1,a2,b2,y0,y1,q){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
     beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*pbeta(q,a1+y1,b1+y0)+
      (1-c)*pbeta(q,a2+y1,b2+y0)
  
  return(result)}

## posterior mean
posterior.mean<-function(a1,b1,a2,b2,y0,y1){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
     beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*(a1+y1)/(a1+b1+y1+y0)+
      (1-c)*(a2+y1)/(a2+b2+y1+y0)
  
  return(result)}

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
    
    futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
    efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)

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
    inner[k,paste0("fut.inf.",time[l])]<-(posterior.quantile(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y0=y0[n[l]],y1=y1[n[l]],q=p.intr)>sig.fut)
    inner[k,paste0("eff.inf.",time[l])]<-(1-posterior.quantile(a1=alpha.skpt,b1=beta.skpt,
                                           a2=alpha.enth,b2=beta.enth,
                                           y0=y0[n[l]],y1=y1[n[l]],q=p.skpt)>sig.eff)
    inner[k,paste0("ss.",time[l])]<-n[l]
    inner[k,paste0("phat.",time[l])]<-y1[n[l]]/n[l]
    inner[k,paste0("post.mean.",time[l])]<-posterior.mean(a1=alpha.skpt,b1=beta.skpt,
                                                     a2=alpha.enth,b2=beta.enth,
                                                     y1=y1[n[l]],y0=y0[n[l]])
    inner[k,paste0("cov.",time[l])]<-posterior.cov(a1=alpha.skpt,b1=beta.skpt,
                                              a2=alpha.enth,b2=beta.enth,
                                              y1=y1[n[l]],y0=y0[n[l]],
                                              p=p.range[j],
                                              cred.tail=cred.tail)
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



par(mfrow = c(1,1)) 
i<-1
plot(p.range,outer[i,,"eff.mon.initial"],type='l',ylim=c(0,1),lwd=2,lty='longdash',
     ylab="Probability",xlab="",
     main="Sequential Design Properties",
     axes=FALSE)
box()
text(p.range,outer[i,,"eff.mon.initial"],labels=format(round(outer[i,,"eff.mon.initial"],digits=2),nsmall=2),pos=3)
lines(p.range,outer[i,,"fut.mon.final"],lwd=2)
lines(p.range,outer[i,,"fut.mon.initial"],lwd=2,lty='longdash')
text(p.range,outer[i,,"fut.mon.initial"],labels=format(round(outer[i,,"fut.mon.initial"],digits=2),nsmall=2),pos=1)
lines(p.range,outer[i,,"eff.mon.final"],lwd=2)


lines(p.range,outer[i,,"inc"],lwd=2,lty='longdash')
lines(p.range,outer[i,,"inc.final"],lwd=2)
#text(p.range,outer[i,,"inc"],
#     labels=format(round(outer[i,,"inc"],digits=2),nsmall=2),pos=1)

axis(1,las=0,at=p.range,labels=format(p.range,nsmall=2))
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')
abline(v=c(p.skpt,p.enth),col='grey',lty='dashed')
row<-1
  for (column in 1:length(p.range)){
    mtext(text=paste0(format(round(outer[k,column,"ss.initial"],digits=1),nsmall=1),
                      " + ",
                      format(round(outer[k,column,"ss.final"]-
                                   outer[k,column,"ss.initial"],digits=1),nsmall=1),
                      " = ",
                      format(round(outer[k,column,"ss.final"],digits=1),nsmall=1)),
                      side=1,line=row+1,at=p.range[column])
  }
row<-2
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",
        format(round(outer[k,column,"post.mean.initial"],digits=3),nsmall=3),
                    " (F) ",
        format(round(outer[k,column,"post.mean.final"],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[column])
}
row<-3
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",
        format(round(outer[k,column,"cov.initial"],digits=3),nsmall=3),
                    " (F) ",
        format(round(outer[k,column,"cov.final"],digits=3),nsmall=3)),
        side=1,line=row+1,at=p.range[column])
}
mtext(text="SS",side=1,line=2,at=0.125)
mtext(text="PM",side=1,line=3,at=0.125)
mtext(text="CP",side=1,line=4,at=0.125)
column<-1
text(locator(1),"Stop Early for Efficacy")
text(locator(1),"Final Efficacy")
text(locator(1),"Stop Early for Futility")
text(locator(1),"Final Futility")
text(locator(1),"Inconclusive")
text(locator(1),"Final Inconclusive")
# legend(x=0.325,y=0.6,
#        legend=c("Stop Early for Efficacy",
#                 "Stop Early for Futility",
#                 "Inconclusive with Full Data"),
#        lty=c('longdash','dotted','solid'),
#        cex=0.5,
#        box.lwd = 1,box.col = "black",bg = "white",pt.cex = 1,
#        text.width=20)

#### END EXAMPLES FOR 7/26/19 ####


## POWER CURVES
par(ask=FALSE)
par(mfrow = c(1,1))
k<-1
plot(p.range,outer.eff[k,],type='l',ylim=c(0,1),lwd=2,
     ylab="Probability",xlab="",
     main="Power Curves (Probability of Efficacy at Interim)\n Modifying Frequency of Monitoring",
     axes=FALSE)
for (k in 2:6){
lines(p.range,outer.eff[k,],lwd=2)
}
box()
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
par(mar=c(5.1+4,4.1,4.1,2.1+1))
par(mfrow = c(1, 1)) 
k1<-4
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff"][out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='4 month follow-up, 2 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
mtext(text="SS",side=1,line=3,at=-.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff"][out.mean==k1 & enr.shape==k2],
     labels=format(round(outer[,,"eff"][out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer[,,"eff.final"][out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer[,,"eff.final"][out.mean==k1 & enr.shape==k2],
     labels=format(round(outer[,,"eff.final"][out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
box()
axis(1,las=0,at=1/freq.mntr[out.mean==k1 & enr.shape==k2],
     labels=80/freq.mntr[out.mean==k1 & enr.shape==k2])
axis(2,las=2,at=seq(0,.15,by=0.01),labels=format(seq(0,.15,by=0.01),nsmall=2))
abline(h=seq(0,0.15,by=0.01),col='grey')
for (column in 1:length(freq.mntr[out.mean==k1 & enr.shape==k2])){
  mtext(text=paste0(format(round(
    outer[,,"ss.initial"][out.mean==k1 & enr.shape==k2][column],
    digits=1),nsmall=1),
    " + ",
    format(round(outer[,,"ss.final"][out.mean==k1 & enr.shape==k2][column]-
                   outer[,,"ss.initial"][out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1),
    " = ",
    format(round(outer[,,"ss.final"][out.mean==k1 & enr.shape==k2][column],
                 digits=1),nsmall=1)),
    side=1,
    line=9-column,
    at=1/freq.mntr[column])
}
k1<-8
k2<-1
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='8 month follow-up, 2 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff.final[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff.final[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
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
    line=9-column,
    at=1/freq.mntr[column])
}
par(mar=c(5.1+4,4.1,4.1,2.1+1))
k1<-4
k2<-0.25
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='4 month follow-up, 8 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff.final[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff.final[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
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
    line=9-column,
    at=1/freq.mntr[column])
}
k1<-8
k2<-0.25
plot(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     type='l',ylim=c(0,.15),lwd=2,lty='longdash',
     ylab="Type 1 Error Rate",xlab="",
     main='8 month follow-up, 8 enrollments per month',#paste0("k1=",k1,", k2=",k2),
     axes=FALSE)
mtext(text="Number of Interim Analyses",
      side=1,
      line=2,
      at=.5)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=3)
lines(1/freq.mntr[out.mean==k1 & enr.shape==k2],
      outer.eff.final[out.mean==k1 & enr.shape==k2],lwd=2)
text(1/freq.mntr[out.mean==k1 & enr.shape==k2],
     outer.eff.final[out.mean==k1 & enr.shape==k2],
     labels=format(round(outer.eff.final[out.mean==k1 & enr.shape==k2],digits=3),
                   nsmall=3),pos=1)
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
    line=9-column,
    at=1/freq.mntr[column])
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
