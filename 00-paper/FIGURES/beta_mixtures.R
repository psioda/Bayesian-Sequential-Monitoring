rm(list = ls())

# mean of skeptical prior
p.skpt<-0.20
# mean of enthuastic prior
p.enth<-0.40
# futility theta
p.intr<-0.30
# value of true response proportion
p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05)
# tail probabilities for priors (low, high)
tail.skpt<-0.045
tail.enth<-0.05
# prior model probabilities
prior.skpt<-1/2 
prior.enth<-1/2 
# maximum sample sizes
max.ss<-80
# significant trial result threshold
sig.fut<-0.85
sig.eff<-0.95
# number of simulated trials per design
reps<-10000
# compute empirical quantities for credible interval
sims<-10000
# credible interval is 1-cred.tail
cred.tail<-0.05

enr.mnths<-c(2)
out.mnths<-c(4)
freq.mntr<-c(38)

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

#################################################################################################
## SIMULATIONS ##################################################################################
#################################################################################################

## Step 1: Create outer loop based on frequency of interim analyses
# stop early for efficacy
outer.eff<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
outer.eff.final<-matrix(nrow=length(freq.mntr),ncol=length(p.range))
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
      print(k) 
      print(i)
      }
    
    efficacy<-0
    futility<-0
    inner.inc[i]<-1
    cutoff.time<-vector()
    event.times<-rep(seq(from=1,
                          to=max.ss/enr.mnths[k],by=1),
                          each=enr.mnths[k])
    responses<-rbinom(n=max.ss,size=1,prob=p.range[j])
    
    for (n in 1:max.ss){
      
      y1<-sum(responses[1:n])
      y0=n-y1
      
      # futility (even optimst would give up)
      futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
      # efficacy (even pessimist would accept)
      efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=FALSE)

      if (n%%freq.mntr[k]==0 & futility>sig.fut){
        inner.fut[i]<-1
        inner.inc[i]<-0
        break
        }
      else if (n%%freq.mntr[k]==0 & efficacy>sig.eff){
        inner.eff[i]<-1
        inner.inc[i]<-0
        break
        }
      else {
        inner.inc[i]<-1
      }
    }
    cutoff.time<-event.times[n]
    responses.final<-responses[event.times<=cutoff.time+out.mnths[k]]
    n.final<-length(responses.final)
    y1.final<-sum(responses.final)
    y0.final<-n.final-y1.final
    
    inner.fut.final[i]<-(pbeta(p.intr,alpha.enth+y1.final,beta.enth+y0.final,
                              lower.tail=TRUE)>sig.fut)
    inner.eff.final[i]<-(pbeta(p.skpt,alpha.skpt+y1.final,beta.skpt+y0.final,
                              lower.tail=FALSE)>sig.eff)
    inner.inc.final[i]<-1-inner.fut.final[i]-inner.eff.final[i]
    
    
    inner.ss.initial[i]<-n
    inner.phat.initial[i]<-y1/n
    inner.ss.final[i]<-n.final
    inner.phat.final[i]<-y1.final/n.final
    
    # inference posterior mixture proportions
    c1.initial<-prior.skpt*beta(alpha.skpt+y1,beta.skpt+y0)/
      (prior.skpt*beta(alpha.skpt+y1,beta.skpt+y0)+
         prior.enth*beta(alpha.enth+y1,beta.enth+y0))
    c2.initial<-prior.enth*beta(alpha.enth+y1,beta.enth+y0)/
      (prior.skpt*beta(alpha.skpt+y1,beta.skpt+y0)+
         prior.enth*beta(alpha.enth+y1,beta.enth+y0))
    c1.final<-prior.skpt*beta(alpha.skpt+y1.final,beta.skpt+y0.final)/
      (prior.skpt*beta(alpha.skpt+y1.final,beta.skpt+y0.final)+
         prior.enth*beta(alpha.enth+y1.final,beta.enth+y0.final))
    c2.final<-prior.enth*beta(alpha.enth+y1.final,beta.enth+y0.final)/
      (prior.skpt*beta(alpha.skpt+y1.final,beta.skpt+y0.final)+
         prior.enth*beta(alpha.enth+y1.final,beta.enth+y0.final))
        # posterior mean
    inner.post.mean.initial[i]<-c1.initial*(alpha.skpt+y1)/(alpha.skpt+beta.skpt+n)+
      c2.initial*(alpha.enth+y1)/(alpha.enth+beta.enth+n)
    inner.post.mean.final[i]<-c1.final*(alpha.skpt+y1.final)/(alpha.skpt+beta.skpt+n.final)+
      c2.final*(alpha.enth+y1.final)/(alpha.enth+beta.enth+n.final)
  
    # theta_1 <- rbeta(round(sims*c1.initial,log10(sims)), alpha.skpt+y1, beta.skpt+y0)
    # theta_2 <- rbeta(round(sims*c2.initial,log10(sims)), alpha.enth+y1, beta.enth+y0)
    # theta <- sort(c(theta_1, theta_2))
    # q_lower <- theta[round((cred.tail / 2) * sims)]
    # q_upper <- theta[round((1 - cred.tail / 2) * sims)]
    # inner.cov.initial[i]=(p.range[j]>q_lower & p.range[j]<q_upper) # initial
    # 
    # theta_1 <- rbeta(round(sims*c1.final,log10(sims)), alpha.skpt+y1, beta.skpt+y0)
    # theta_2 <- rbeta(round(sims*c2.final,log10(sims)), alpha.enth+y1, beta.enth+y0)
    # theta <- sort(c(theta_1, theta_2))
    # q_lower <- theta[round((cred.tail / 2) * sims)]
    # q_upper <- theta[round((1 - cred.tail / 2) * sims)]
    # inner.cov.final[i]=(p.range[j]>q_lower & p.range[j]<q_upper) # final
    
  }

  outer.ss.initial[k,j]<-mean(inner.ss.initial)
  outer.ss.final[k,j]<-mean(inner.ss.final)
  outer.phat.initial[k,j]<-mean(inner.phat.initial)
  outer.phat.final[k,j]<-mean(inner.phat.final)
  outer.fut[k,j]<-mean(inner.fut)
  outer.fut.final[k,j]<-mean(inner.fut.final)
  outer.eff[k,j]<-mean(inner.eff)
  outer.eff.final[k,j]<-mean(inner.eff.final)
  outer.inc[k,j]<-mean(inner.inc)
  outer.inc.final[k,j]<-mean(inner.inc.final)
  outer.post.mean.initial[k,j]<-mean(inner.post.mean.initial)
  outer.post.mean.final[k,j]<-mean(inner.post.mean.final)
  outer.cov.initial[k,j]<-mean(inner.cov.initial)
  outer.cov.final[k,j]<-mean(inner.cov.final)
}
}
outer.ss.initial-outer.ss.final
outer.cov.initial-outer.cov.final

k<-1
plot(p.range,outer.fut[k,],type='l',ylim=c(0,1),col='red',lwd=2,lty='longdash',
     ylab="Probability",xlab="",main="Sequential Design Properties",
     axes=FALSE)
box()
lines(p.range,outer.fut.final[k,],col='red',lwd=2)
lines(p.range,outer.eff[k,],lwd=2,lty='longdash',col='green')
lines(p.range,outer.eff.final[k,],lwd=2,col='green')
lines(p.range,outer.inc[k,],lwd=2,lty='longdash')
lines(p.range,outer.inc.final[k,],lwd=2)
axis(1,las=0,at=p.range,labels=format(p.range,nsmall=2))
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')
abline(v=c(p.skpt,p.enth),col='grey',lty='dashed')
row<-1
  for (column in 1:length(p.range)){
    mtext(text=paste0(round(outer.ss.initial[k,column],digits=1),
                      " + ",round(outer.ss.final[k,column]-
                                  outer.ss.initial[k,column],digits=1),
                      " = ",round(outer.ss.final[k,column],digits=1)),
                      side=1,line=row+1,at=p.range[column])
  }
row<-2
for (column in 1:length(p.range)){
  mtext(text=paste0("(I) ",round(outer.post.mean.initial[k,column],digits=3),
                    " (F) ",round(outer.post.mean.final[k,column],digits=3)),
        side=1,line=row+1,at=p.range[column])
}
column<-1
legend(x=0.35,y=0.5,legend=c("Stop Early for Efficacy","Stop Early for Futility",
                "Inconclusive w/ Full Dataset"),
       col=c("green","red","black"),lty=c('longdash','dotted','solid'),
       cex=0.8,
       box.lwd = 1,box.col = "black",bg = "white",pt.cex = 1)
mtext(text="SS",side=1,line=2,at=0.125)
mtext(text="PM",side=1,line=3,at=0.125)
#### END EXAMPLES FOR 7/26/19 ####

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
plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",col="red",xlab="Response Probability",ylab="Density Value",main="Skeptical Prior",
     #xaxt="n",
     ylim=c(0,max(dbeta(x,alpha.skpt,beta.skpt),dbeta(x,alpha.enth,beta.enth))))
#axis(1,at=p.skpt,labels=expression(theta[0]))
#axis(1,at=seq(0,1,by=0.2))
abline(v=p.enth)
abline(v=p.skpt)
#high (enthuastic)
plot(x,dbeta(x,alpha.enth,beta.enth),type="l",col="blue",
     xlab="Response Probability",ylab="Density Value",main="Enthuastic Prior",
     #xaxt="n",
     ylim=c(0,max(dbeta(x,alpha.skpt,beta.skpt),dbeta(x,alpha.enth,beta.enth))))
#axis(1,at=p.enth,labels=expression(theta[A]))
abline(v=p.enth)
abline(v=p.skpt)
# Step 6: Plot inference priors (mixtures)
## 7/19/19 make 50:50 solid black line, 25:75 dark grey with different dash patterns,
# make skeptical/enthuastic light grey with smaller thinkness
# no legend, put omega on actual lines

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
