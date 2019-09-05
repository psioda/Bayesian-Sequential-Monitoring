## The purpose of this program is to find parameterize 
## mixtures of beta distributions for spike/slab and 
## flat versions of skeptical and enthuastic priors

require(flexmix)
require(betareg)

invcdf.lower<-function(x,t){t*sqrt(x)}
invcdf.upper<-function(x,t){1-(1-t)*sqrt(1-x)}
invcdf.lower.scaled<-function(x,t){t*sqrt(x/t)}
invcdf.upper.scaled<-function(x,t){1-(1-t)*sqrt((1-x)/(1-t))}
#plot(x,invcdf.lower.scaled(x,t=t.lower))
#plot(x,invcdf.upper.scaled(x,t=t.upper))
#hist(invcdf.lower.scaled(x,t=t.lower))
#hist(invcdf.upper.scaled(x,t=t.upper))

## PARAMETERS ##
t.lower<-0.2
target.cutoff<-t.lower 
target.mean<-0.4
t.upper<-0.7
cut<-c(0.3,0.5)        # cutoff points determining inner 3 segments
target.tail<-0.05        # created with lower tail in mind
p<-c(target.tail,0.225,NA,NA,0.01)  # probabilities of the 5 segments

mean.tol<-0.01           
tail.tol<-0.01

E<-c(t.lower*(2/3),       # expected value of each segment separately
     (t.lower+cut[1])/2,  # used the fact that the end pieces are triangular
     (cut[1]+cut[2])/2,
     (cut[2]+t.upper)/2,
     t.upper+(1-t.upper)*1/3)

# solve linear equation
A<-t(matrix(data=c(1,1,(cut[1]+cut[2])/2,(cut[2]+t.upper)/2),
          nrow=2,ncol=2))
b<-c(1-sum(p,na.rm=TRUE),target.mean-sum(E*p,na.rm=TRUE))
p[3]<-solve(a=A,b=b)[1]
p[4]<-solve(a=A,b=b)[2]

# reps=1E3
# samps<-c(invcdf.lower.scaled(x=seq(0,t.lower,length=reps*p[1]),t=t.lower),
#          seq(t.lower,cut[1],length=reps*p[2]),
#          seq(cut[1],cut[2],length=reps*p[3]),
#          seq(cut[2],t.upper,length=reps*p[4]),
#          invcdf.upper.scaled(x=seq(t.upper,1,length=reps*p[5]),t=t.upper))

# check histogram
reps<-1E6
samps<-c(invcdf.lower(x=runif(n=reps*p[1]),t=t.lower),
         t.lower+runif(n=reps*p[2])*(cut[1]-t.lower),
         cut[1]+runif(n=reps*p[3])*(cut[2]-cut[1]),
         cut[2]+runif(n=reps*p[4])*(t.upper-cut[2]),
         invcdf.upper(x=runif(n=reps*p[5]),t=t.upper))
hist(samps,breaks=seq(0,1,length=100),freq=F)
mean(samps)

cond1<-1
cond2<-1

outer<- array(NA,dim = c(1000,4,3),
              dimnames = list(seq_len(1000),seq_len(4),c("w","a","b")))


#while(cond1>mean.tol | cond2 > tail.tol){
for (i in 1:1){ 
print(i)
  
reps<-1000
# samps<-c(invcdf.lower(x=runif(n=reps*p[1]),t=t.lower),
#          t.lower+runif(n=reps*p[2])*(cut[1]-t.lower),
#          cut[1]+runif(n=reps*p[3])*(cut[2]-cut[1]),
#          cut[2]+runif(n=reps*p[4])*(t.upper-cut[2]),
#          invcdf.upper(x=runif(n=reps*p[5]),t=t.upper))
samps<-c(invcdf.lower.scaled(x=seq(0,t.lower,length=reps*p[1]),t=t.lower),
         seq(t.lower,cut[1],length=reps*p[2]),
         seq(cut[1],cut[2],length=reps*p[3]),
         seq(cut[2],t.upper,length=reps*p[4]),
         invcdf.upper.scaled(x=seq(t.upper,1,length=reps*p[5]),t=t.upper))
hist(samps,breaks=seq(0,1,length=100),freq=F)
mean(samps)

d <- data.frame(y=samps[-c(1,length(samps))])

k<-4
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

index<-1
x<-seq(0,1,by=0.01)
hist(samps,breaks=seq(0,1,length=100),freq=F)
plot(xlim=c(0,1),ylim=c(0,1))
par(mfrow=c(1,1))

lines(x,outer[index,1,"w"]*dbeta(x,shape1=outer[index,1,"a"], shape2=outer[index,1,"b"])+
       outer[index,2,"w"]*dbeta(x,shape1=outer[index,2,"a"], shape2=outer[index,2,"b"])+
       outer[index,3,"w"]*dbeta(x,shape1=outer[index,3,"a"], shape2=outer[index,3,"b"])+
       outer[index,4,"w"]*dbeta(x,shape1=outer[index,4,"a"], shape2=outer[index,4,"b"]))
lines(x,dbeta(x,alpha.enth,beta.enth))


outer[,,"w"]*outer[,,"a"]/(outer[,,"a"]/outer[,,"b"])
sum(w*a/(a+b))
sum(w*pbeta(q=target.cutoff,shape1=a,shape2=b))
x<-seq(0,1,by=0.01)
hist(d$y,freq=FALSE,ylim=c(0,10))
lines(x, dbeta(x, shape1 = a[1], shape2 = b[1]),col = hcl(0, 80, 50), lwd = 2)
lines(x, dbeta(x, shape1 = a[2], shape2 = b[2]),col = hcl(240, 80, 50), lwd = 2)
lines(x, w[1]*dbeta(x, shape1 = a[1], shape2 = b[1])+w[2]*dbeta(x, shape1 = a[2], shape2 = b[2]))
lines(x,dbeta(x,alpha.skpt,beta.skpt))

test1<-rowSums(outer[,,"w"]*outer[,,"a"]/(outer[,,"a"]+outer[,,"b"]))
test2<-rowSums(outer[,,"w"]*pbeta(q=target.cutoff,
                                  shape1=outer[,,"a"],
                                  shape2=outer[,,"b"]))
test3<-test2-0.05
head(sort(test3),n=30)
index<-126 
plot(x, outer[index,1,"w"]*dbeta(x,shape1=outer[index,1,"a"], shape2=outer[index,1,"b"])+
       outer[index,2,"w"]*dbeta(x,shape1=outer[index,2,"a"], shape2=outer[index,2,"b"]))
x<-seq(0,1,by=0.01)
plot(x, outer[index,1,"w"]*dbeta(x,shape1=outer[index,1,"a"], shape2=outer[index,1,"b"])+
        outer[index,2,"w"]*dbeta(x,shape1=outer[index,2,"a"], shape2=outer[index,2,"b"]))

test2<-
test2<-rowSums(test1)
mean(test2)
