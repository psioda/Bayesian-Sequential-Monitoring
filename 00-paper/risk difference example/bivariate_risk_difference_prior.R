rm(list = ls())
require(rmutil)
require(lattice)
require(pracma)

sigma<-4
alpha<-2
delta.enth<-0.2
delta.skpt<-0

# graph bivariate skeptical prior
denB<-function(x,y){
  exp(-((sigma*((x-y)-delta.skpt))^2)^alpha)
}
denB.nc<-int2(denB,a=c(0,0),b=c(1,1))
x <- seq(0, 1, length= 100)
y <- x
z <- outer(x, y, denB)
wireframe(z, drape=T, col.regions=rainbow(100))

# sample data
y1.x<-c(0,5,10,20)
y0.x<-c(0,5,10,20)
y1.y<-c(0,5,10,20)
y0.y<-c(0,5,10,20)

results<-NA
exp.x<-NA
exp.y<-NA
greater<-NA

for (i in 1:length(y1.x)){
  
posterior.skpt<-function(x,y){
  exp(y1.x[i]*log(x)+y0.x[i]*log(1-x))*
  exp(y1.y[i]*log(y)+y0.y[i]*log(1-y))*
  exp(-((sigma*((x-y)-delta.skpt))^2)^alpha)
}
posterior.skpt.nc<-integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1)[[1]]

# check
results[i]<-integral2(posterior.skpt,xmin=0,xmax=1,ymin=0,ymax=1)[[1]]/posterior.skpt.nc

# skeptical expectation x
posterior.x<-function(x,y){
  x*
  exp(y1.x[i]*log(x)+y0.x[i]*log(1-x))*
  exp(y1.y[i]*log(y)+y0.y[i]*log(1-y))*
    exp(-((sigma*((x-y)-delta.skpt))^2)^alpha)
}
# skeptical expectation y
posterior.y<-function(x,y){
  y*
  exp(y1.x[i]*log(x)+y0.x[i]*log(1-x))*
  exp(y1.y[i]*log(y)+y0.y[i]*log(1-y))*
  exp(-((sigma*((x-y)-delta.skpt))^2)^alpha)
}

exp.x[i]<-integral2(posterior.x,xmin=0,xmax=1,ymin=0,ymax=1)[[1]]/posterior.skpt.nc
exp.y[i]<-integral2(posterior.y,xmin=0,xmax=1,ymin=0,ymax=1)[[1]]/posterior.skpt.nc

# integrate over different areas
ymax <- function(x) x-0.1
greater[i]<-integral2(posterior.skpt,xmin=0.1,xmax=1,ymin=0,ymax)[[1]]/posterior.skpt.nc
}

## SIMULATIONS ##################################################################################


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


for (i in 1:length(freq.mntr)){
  for (j in 1:length(p.range)){
    for (k in 1:reps){
      
      if (k%%1000==0){print(paste0("Model ",i,", Response p=",p.range[j],", Simulation ",k))}
      
      responses.x<-rep(c(0,1),max.ss/2)
      responses.y<-responses.x
      y1.x<-cumsum(responses.x)
      y1.y<-cumsum(responses.y)
      y0.x<-seq(1:length(responses.x))-y1.x
      y0.y<-seq(1:length(responses.y))-y1.y
      
      futility<-rep(NA,length=max.ss)
      efficacy<-rep(NA,length=max.ss)
      
      for (h in 1:max.ss){
        # skeptical posterior & efficacy
        posterior.skpt<-function(x,y){exp(y1.x*log(x)+y0.x*log(1-x))*
                                      exp(y1.y*log(y)+y0.y*log(1-y))*
                                      exp(-((sigma*((x-y)-delta.skpt))^2)^alpha)}
        posterior.skpt.nc<-int2(posterior.skpt,a=c(0,0),b=c(1,1))[[1]]
        
        posterior.x<-function(x,y){
          x*
            exp(y1.x*log(x)+y0.x*log(1-x))*
            exp(y1.y*log(y)+y0.y*log(1-y))*
            exp(-((sigma*((x-y)-delta))^2)^alpha)
        }
        posterior.y<-function(x,y){
          y*
            exp(y1.x*log(x)+y0.x*log(1-x))*
            exp(y1.y*log(y)+y0.y*log(1-y))*
            exp(-((sigma*((x-y)-delta))^2)^alpha)
        }
        
        int2(posterior.x,a=c(0,0),b=c(1,1))/posterior.nc[[1]]
        int2(posterior.y,a=c(0,0),b=c(1,1))/posterior.nc[[1]]
        
          posterior.skpt<-function(x){
            exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
          }
          nc.skpt<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
          posterior.nc.skpt<-function(x){
            exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/nc.skpt
          }
          efficacy[h]<-integrate(posterior.nc.skpt,lower=0+epsilon,upper=p.skpt)[[1]]
          
          # enthuastic posterior & futility
          posterior.enth<-function(x){
            exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
          }
          nc.enth<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
          posterior.nc.enth<-function(x){
            exp(y1[h]*log(x)+y0[h]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/nc.enth
          }
          futility[h]<-integrate(posterior.nc.enth,lower=0+epsilon,upper=p.intr)[[1]]
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
    }
  }
}


posterior.x<-function(x,y){
  x*
    exp(y1.x*log(x)+y0.x*log(1-x))*
    exp(y1.y*log(y)+y0.y*log(1-y))*
    exp(-((sigma*((x-y)-delta))^2)^alpha)
}
posterior.y<-function(x,y){
  y*
    exp(y1.x*log(x)+y0.x*log(1-x))*
    exp(y1.y*log(y)+y0.y*log(1-y))*
    exp(-((sigma*((x-y)-delta))^2)^alpha)
}



##  Compute the volume of a sphere
f <- function(x, y) sqrt(1 -x^2 - y^2)
xmin <- 0; xmax <- 1
ymin <- 0; ymax <- function(x) sqrt(1 - x^2)
integral2(f, xmin, xmax, ymin, ymax)
##  Integrate 1/( sqrt(x + y)*(1 + x + y)^2 ) over the triangle
##   0 <= x <= 1, 0 <= y <= 1 - x.  The integrand is infinite at (0,0).
f <- function(x,y) 1/( sqrt(x + y) * (1 + x + y)^2 )
ymax <- function(x) 1 - x
integral2(f, 0,1, 0,ymax)

f <-function(x,y) 1
ymax<-function(x) x
integral2(f,0,1,0,ymax)
