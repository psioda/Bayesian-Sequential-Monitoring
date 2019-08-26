rm(list = ls())
require(flexmix)
require(betareg)

#################################################################################################
## Generate fake data ###########################################################################
#################################################################################################

# # parameters of distribution #1
# alpha1 <- 10
# beta1 <- 30
# # parameters of distribution #2
# alpha2 <- 30
# beta2 <- 10
# # Generate bimodal data
# set.seed(0)
# d <- data.frame(y = c(rbeta(100, alpha1, beta1), rbeta(50, alpha2, beta2)))
# n <- length(d$y)
# # histogram
# par(mfrow=c(1,1))
# hist(d$y, 50)

# q<-c(0,0.15,0.25,0.40,1.00)
# p<-c(0,0.205,0.70,0.045,0.05)
# 
# # do a grid search for this piece
# 
# (sum(p)==1) # check p
# x<-seq(0,1,length=10000)
# x<-x[-c(1,length(x))]
# fx<-0
# # create the distribution function
# for (i in 1:length(x)){
#   if (q[1]<x[i] & x[i]<=q[2]) fx[i]=p[2]/(q[2]-q[1])
#   if (q[2]<x[i] & x[i]<=q[3]) fx[i]=p[3]/(q[3]-q[2])
#   if (q[3]<x[i] & x[i]<=q[4]) fx[i]=p[4]/(q[4]-q[3])
#   if (q[4]<x[i] & x[i]<=q[5]) fx[i]=2*p[5]*(1-x[i])/(q[5]-q[4])^2 # make last segment a triangle
# }
# make data to match the distribution function

# # sampling from the lower triangular section is a little tricky
# # https://en.wikipedia.org/wiki/Inverse_transform_sampling
# x <- seq(0,1,by=0.1)
# plot(x,2*x,type='l') # pdf 
# plot(x,x^2,type='l') # cdf
#  invcdf<-function(x){
#   sqrt(x)
# }
# plot(x,invcdf(x),type='l') # inverse cdf
# hist(invcdf(runif(100000,0,1)),freq=F) # random draws

# x<-seq(0.4,1,by=0.1)
# plot(x, (1-x)/0.6^2*2,type='l') # pdf
# cdf<-function(x){
#   -1.77778 + 5.55556*x - 2.77778*x^2
# }
# plot(x,cdf(x),type='l') # cdf
# invcdf<-function(x){
#   1 - 0.6*sqrt(1 - x)
# }
# plot(x,invcdf(x),type='l') # inverse cdf
# U<-runif(100000,0,1)
# hist(invcdf(U),freq=F) # random draws



# CDF for lower triangular portion
invcdf<-function(x){
  1 - 0.6*sqrt(1 - x)
}
reps<-1E6
samps<-c(seq(0,0.15,length=reps*0.205),
  seq(0.15,0.25,length=reps*0.7),
  seq(0.25,0.40,length=reps*0.045),
  invcdf(seq(0,1,length=reps*0.05)))
hist(samps)

d <- data.frame(y=samps[-c(1,length(samps))])

#################################################################################################
## Gitting mixtures of beta distributions ######################################################
#################################################################################################

m <- betamix(y ~ 1 | 1, data = d, k = 2)
mu <- plogis(coef(m)[,1]) # plogis distribution of logistic function # [,1] is intercept
phi <- exp(coef(m)[,2])   # [,2] is phi
a <- mu * phi
b <- (1 - mu) * phi 
w <- prior(m$flexmix) # probability of clusters for mixture


sum(w*a/(a+b)) # find posterior mean of mixture
mean(d$y) # compare with mean of data
sum(w*pbeta(q=0.4,shape1=a,shape2=b)) # find desired tail probability of mixture

x<-seq(0,1,by=0.01)
## lines for fitted densities
hist(d$y,freq=FALSE)
lines(x, dbeta(x, shape1 = a[1], shape2 = b[1]),col = hcl(0, 80, 50), lwd = 2)
lines(x, dbeta(x, shape1 = a[2], shape2 = b[2]),col = hcl(240, 80, 50), lwd = 2)
lines(x, w[1]*dbeta(x, shape1 = a[1], shape2 = b[1])+w[2]*dbeta(x, shape1 = a[2], shape2 = b[2]))

