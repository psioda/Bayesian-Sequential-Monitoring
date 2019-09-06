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

inner <- array(NA, dim=c(reps,length(freq.mntr),length(matrix.names)), 
               dimnames = list(seq_len(reps),freq.mntr,matrix.names))

for (j in 1:length(p.range)){
  
  for (k in 1:reps){
    
    if (k%%1000==0){print(paste0("Response p=",p.range[j],", Simulation ",k))}
    
enr.times<-cumsum(rgamma(n=max.ss,shape=enr.shape,scale=0.5))
outcome.times<-sort(enr.times+rnorm(n=max.ss,mean=out.mean,sd=0.25))
responses<-rbinom(n=max.ss,size=1,prob=p.range[j])
y1<-cumsum(responses)
y0<-seq(1:length(responses))-y1
futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)

for (i in 1:length(freq.mntr)){
   		n.initial<-min(
   		 which(((futility>sig.fut) | ((1-efficacy)>sig.eff)) & (seq(1:max.ss)%%freq.mntr[i]==0)),
   		max.ss,na.rm=TRUE)
cutoff.time<-outcome.times[n.initial]
n.final<-length(responses[enr.times<=cutoff.time])
    		n<-c(n.initial,n.final)

time<-c("initial","final")
    			for (l in 1:2){
   			inner[k,i,paste0("fut.mon.",time[l])]<-(futility[n[l]]>sig.fut)
   			inner[k,i,paste0("eff.mon.",time[l])]<-(1-efficacy[n[l]]>sig.eff)
inner[k,i,paste0("ss.",time[l])]<-n[l]
inner[k,i,paste0("phat.",time[l])]<-y1[n[l]]/n[l]
inner[k,i,paste0("fut.inf.",time[l])]<-(posterior.cdf(
a1=alpha.skpt,b1=beta.skpt,a2=alpha.enth,b2=beta.enth,
y0=y0[n[l]],y1=y1[n[l]],q=p.intr,w=0.5)>sig.fut)
inner[k,i,paste0("eff.inf.",time[l])]<-(1-posterior.cdf(
a1=alpha.skpt,b1=beta.skpt,a2=alpha.enth,b2=beta.enth,
y0=y0[n[l]],y1=y1[n[l]],q=p.skpt,w=0.5)>sig.eff)
    			inner[k,i,paste0("post.mean.",time[l])]<-posterior.mean(
a1=alpha.skpt,b1=beta.skpt,a2=alpha.enth,b2=beta.enth,
y1=y1[n[l]],y0=y0[n[l]])
inner[k,i,paste0("cov.",time[l])]<-posterior.cov(
a1=alpha.skpt,b1=beta.skpt,a2=alpha.enth,b2=beta.enth,
y1=y1[n[l]],y0=y0[n[l]],p=p.range[j], cred.tail=cred.tail)
} # end L LABELS
} # end I FREQ.MNTR
	} # end K REPS
outer[,j,]<-apply(inner,MARGIN=c(2,3),FUN=mean)
} # end J P.RANGE
