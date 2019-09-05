## The purpose of this program is to find parameterize 
## mixtures of beta distributions for spike/slab and 
## flat versions of skeptical and enthuastic priors

require(flexmix)
require(betareg)

invcdf.lower<-function(x,t){t*sqrt(x)}
invcdf.upper<-function(x,t){1-(1-t)*sqrt(1-x)}

## PARAMETERS ##
t.lower<-0.2
target.cutoff<-t.lower 
target.mean<-0.4
t.upper<-0.5625
cut<-0.5
target.tail<-0.05        # created with lower tail in mind
p<-c(target.tail,NA,NA,0.15)  # probabilities of the 5 segments

mean.tol<-0.01           
tail.tol<-0.01

E<-c(t.lower*(2/3),       # expected value of each segment separately
     (t.lower+cut)/2,  # used the fact that the end pieces are triangular
     (cut+t.upper)/2,
     t.upper+(1-t.upper)*1/3)

# solve linear equation
A<-t(matrix(data=c(1,1,(t.lower+cut)/2,(cut+t.upper)/2),
            nrow=2,ncol=2))
b<-c(1-sum(p,na.rm=TRUE),target.mean-sum(E*p,na.rm=TRUE))
p[2]<-solve(a=A,b=b)[1]
p[3]<-solve(a=A,b=b)[2]

# check histogram
reps<-1E6
samps<-c(invcdf.lower(x=runif(n=reps*p[1]),t=t.lower),
         t.lower+runif(n=reps*p[2])*(cut-t.lower),
         cut+runif(n=reps*p[3])*(t.upper-cut),
         invcdf.upper(x=runif(n=reps*p[4]),t=t.upper))
hist(samps,breaks=seq(0,1,length=100),freq=F)
mean(samps)


outer<- array(NA,dim = c(1000,3,3),
              dimnames = list(seq_len(1000),seq_len(3),c("w","a","b")))
#while(cond1>mean.tol | cond2 > tail.tol){
for (i in 1:1000){ 

  print(i)
  
  reps<-750
  samps<-c(invcdf.lower(x=runif(n=reps*p[1]),t=t.lower),
           t.lower+runif(n=reps*p[2])*(cut-t.lower),
           cut+runif(n=reps*p[3])*(t.upper-cut),
           invcdf.upper(x=runif(n=reps*p[4]),t=t.upper))
  hist(samps,breaks=seq(0,1,length=100),freq=F)
  mean(samps)
  
  d <- data.frame(y=samps[-c(1,length(samps))])
  
  k<-3
  m <- betamix(y ~ 1 | 1, data = d, k = k)
  
  if (length(coef(m))==k*2){
    mu <- plogis(coef(m)[,1]) # plogis distribution of logistic function # [,1] is intercept
    phi <- exp(coef(m)[,2])   # [,2] is phi
    a <- mu * phi
    b <- (1 - mu) * phi 
    w <- prior(m$flexmix) # probability of clusters for mixture
    cond1<-abs(target.mean-sum(w*a/(a+b))) # find posterior mean of mixture
    cond2<-abs(target.tail-sum(w*pbeta(q=target.cutoff,shape1=a,shape2=b))) # find desired tail probability of mixture
    print(cond1)
    print(cond2)
    outer[i,,"w"]<-w
    outer[i,,"a"]<-a
    outer[i,,"b"]<-b
    }
}


index<-40 
x<-seq(0,1,by=0.01)
plot(x, outer[index,1,"w"]*dbeta(x,shape1=outer[index,1,"a"], shape2=outer[index,1,"b"])+
        outer[index,2,"w"]*dbeta(x,shape1=outer[index,2,"a"], shape2=outer[index,2,"b"])+
        outer[index,3,"w"]*dbeta(x,shape1=outer[index,3,"a"], shape2=outer[index,3,"b"]),
     ylim=c(0,2.5))

sum(w*a/(a+b))
sum(w*pbeta(q=target.cutoff,shape1=a,shape2=b))
x<-seq(0,1,by=0.01)
hist(d$y,freq=FALSE,ylim=c(0,10),xlim=c(0,1))
#lines(x, dbeta(x, shape1 = a[1], shape2 = b[1]),col = hcl(0, 80, 50), lwd = 2)
#lines(x, dbeta(x, shape1 = a[2], shape2 = b[2]),col = hcl(240, 80, 50), lwd = 2)
#lines(x, dbeta(x, shape1 = a[3], shape2 = b[3]),col = hcl(240, 80, 50), lwd = 2)
lines(x, w[1]*dbeta(x, shape1 = a[1], shape2 = b[1])+
         w[2]*dbeta(x, shape1 = a[2], shape2 = b[2])+
         w[3]*dbeta(x, shape1 = a[3], shape2 = b[3])+
         w[4]*dbeta(x, shape1 = a[4], shape2 = b[4])
      )


#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################
p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
p.enth<-0.40
p.intr<-0.30
tail.skpt<-0.045  # tail probabilities for priors (low, high)
tail.enth<-0.05

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

pbeta(q=target.cutoff,shape1=alpha.enth,shape2=beta.enth)
lines(x,dbeta(x,alpha.enth,beta.enth))

