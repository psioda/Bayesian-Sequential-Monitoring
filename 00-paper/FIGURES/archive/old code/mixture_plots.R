## NOT NEEDED ANYMORE
## ALL RELEVANT CODE WENT INTO spike_slab_parameterization

rm(list = ls())

### ORIGINAL FIRST
# mean of skeptical prior
p.skpt<-0.20
# mean of enthuastic prior
p.enth<-0.40
tail.skpt<-0.045
tail.enth<-0.05

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,1000,by=0.01)

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



# Step 5: Plot Skeptical and enthuastic priors separately
## 7/19/19 match slide 14 from FDA presentation, shade in gray scale, set 350 DPI
par(ask=TRUE)
par(mfrow = c(1, 1)) 
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

### SPIKE NEXT
# mean of skeptical prior
p.skpt<-0.20
# mean of enthuastic prior
p.enth<-0.40
tail.skpt<-0.045*2
tail.enth<-0.05*2

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,1000,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.skpt
quantiles.skpt<-qbeta(tail.skpt,(p.skpt)*phi.range,(1-(p.skpt))*phi.range,lower.tail=FALSE)
# lower tail probability equal to tail.enth
quantiles.enth<-qbeta(tail.enth,(p.enth)*phi.range,(1-(p.enth))*phi.range,lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi.range[which.min(abs(p.enth-quantiles.skpt))] # fixed 5/13/19
phi_H<-phi.range[which.min(abs(p.skpt-quantiles.enth))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.skpt.1<-(p.skpt)*phi_L
beta.skpt.1<-(1-(p.skpt))*phi_L
alpha.enth.1<-(p.enth)*phi_H
beta.enth.1<-(1-(p.enth))*phi_H

par(ask=FALSE)
par(mfrow = c(1, 1)) 
x<-seq(0,1,by=0.01)
const<-20


y<-.5*dbeta(x,alpha.skpt.1,beta.skpt.1)+.5*dbeta(x,1*const,4*const)
plot(x,y,type='l',
     xlab="Response Probability",ylab="Density Value",
     main="Spike-Slab Skeptical Beta Prior")
polygon(c(x[x<=0.4],0.4),c(y[x<=0.4],0),col="blue")
polygon(c(x[x>=0.4],0.4),c(y[x>=0.4],0),col="red")
lines(x,dbeta(x,alpha.skpt,beta.skpt),type='l',lty="dashed",lwd=2)

.5*pbeta(.4,alpha.skpt.1,beta.skpt.1)+.5*pbeta(.4,1*const,4*const)

plot(x,.5*dbeta(x,alpha.enth.1,beta.enth.1)+.5*dbeta(x,40,60),type='l',lwd=2,
     xlab="Response Probability",ylab="Density Value",
     main="Spike-Slab Enthuastic Beta Prior")
y<-.5*dbeta(x,alpha.enth.1,beta.enth.1)+.5*dbeta(x,40,60)
polygon(c(x[x<=0.2],0.2),c(y[x<=0.2],0),col="red")
polygon(c(x[x>=0.2],0.2),c(y[x>=0.2],0),col="blue")
lines(x,dbeta(x,alpha.enth,beta.enth),type='l',lty="dashed",lwd=2)
.5*pbeta(.2,alpha.enth.1,beta.enth.1)+.5*pbeta(.2,40,60)
### FLAT LAST
# mean of skeptical prior
separation<-0.06
# specify the ratio for how much lower vs. higher
# 2:1
p.skpt.1<-0.2-separation
p.skpt.2<-0.2+separation

tail.skpt<-0.045

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,1000,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.skpt
quantiles.skpt.1<-qbeta(tail.skpt-0.03,(p.skpt.1)*phi.range,
                        (1-(p.skpt.1))*phi.range,lower.tail=FALSE)
quantiles.skpt.2<-qbeta(tail.skpt+0.03,(p.skpt.2)*phi.range,
                        (1-(p.skpt.2))*phi.range,lower.tail=FALSE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L.1<-phi.range[which.min(abs(p.enth-quantiles.skpt.1))] # fixed 5/13/19
phi_L.2<-phi.range[which.min(abs(p.enth-quantiles.skpt.2))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.skpt2.1<-(p.skpt.1)*phi_L.1
beta.skpt2.1<-(1-(p.skpt.1))*phi_L.1
alpha.skpt2.2<-(p.skpt.2)*phi_L.2
beta.skpt2.2<-(1-(p.skpt.2))*phi_L.2

y<-1/2*dbeta(x,alpha.skpt2.1,beta.skpt2.1)+
  1/2*dbeta(x,alpha.skpt2.2,beta.skpt2.2)
plot(x,dbeta(x,alpha.skpt,beta.skpt),type='l',
     xlab="Response Probability",ylab="Density Value",
     main="Dampened Skeptical Beta Prior",lty="dashed",lwd=2)
lines(x,y)
polygon(c(x[x<=0.4],0.4),c(y[x<=0.4],0),col="blue")
polygon(c(x[x>=0.4],0.4),c(y[x>=0.4],0),col="red")
lines(x,dbeta(x,alpha.skpt,beta.skpt),type='l',lty="dashed",lwd=2)

pbeta(0.4,alpha.skpt,beta.skpt,lower.tail=FALSE)
1/2*pbeta(0.4,alpha.skpt2.1,beta.skpt2.1,lower.tail=FALSE)+
  1/2*pbeta(0.4,alpha.skpt2.2,beta.skpt2.2,lower.tail=FALSE)

### FLAT LAST
# mean of enthuastic prior
separation<-0.09
p.enth.1<-0.4-separation
p.enth.2<-0.4+separation

tail.enth<-0.05

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,1000,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.enth
quantiles.enth.1<-qbeta(tail.enth+0.049,(p.enth.1)*phi.range,
                        (1-(p.enth.1))*phi.range,lower.tail=TRUE)
quantiles.enth.2<-qbeta(tail.enth-0.049,(p.enth.2)*phi.range,
                        (1-(p.enth.2))*phi.range,lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_H.1<-phi.range[which.min(abs(p.skpt-quantiles.enth.1))] # fixed 5/13/19
phi_H.2<-phi.range[which.min(abs(p.skpt-quantiles.enth.2))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.enth2.1<-(p.enth.1)*phi_H.1
beta.enth2.1<-(1-(p.enth.1))*phi_H.1
alpha.enth2.2<-(p.enth.2)*phi_H.2
beta.enth2.2<-(1-(p.enth.2))*phi_H.2

plot(x,dbeta(x,alpha.enth,beta.enth),type='l',lty="dashed",lwd=2,
     xlab="Response Probability",ylab="Density Value",
     main="Dampened Enthuastic Beta Prior")
y<-1/2*dbeta(x,alpha.enth2.1,beta.enth2.1)+
  1/2*dbeta(x,alpha.enth2.2,beta.enth2.2)
lines(x,y)
polygon(c(x[x<=0.2],0.2),c(y[x<=0.2],0),col="red")
polygon(c(x[x>=0.2],0.2),c(y[x>=0.2],0),col="blue")
lines(x,dbeta(x,alpha.enth,beta.enth),type='l',lty="dashed",lwd=2)





pbeta(0.2,alpha.enth,beta.enth,lower.tail=TRUE)
1/2*pbeta(0.2,alpha.enth2.1,beta.enth2.1,lower.tail=TRUE)+
  1/2*pbeta(0.2,alpha.enth2.2,beta.enth2.2,lower.tail=TRUE)

.5*alpha.enth2.1/(alpha.enth2.1+beta.enth2.1)+
  .5*alpha.enth2.2/(alpha.enth2.2+beta.enth2.2)


