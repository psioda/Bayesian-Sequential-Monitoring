efficacy<-NA
futility<-NA
mean<-NA
lower<-NA
upper<-NA
for (i in 1:length(n)){
  efficacy[i]<-eff_fut(i)[1]
  futility[i]<-eff_fut(i)[2]
  mean[i]<-pm_cp(i)[1]
  lower[i]<-pm_cp(i)[4]
  upper[i]<-pm_cp(i)[5]
}

dens<-array(NA,dim=c(length(n),2,length(grid)),
            dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
loc<-array(NA,dim=c(length(n),2,length(grid)),
           dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
x<-array(NA,dim=c(length(n),2,length(grid)),
         dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
m<-array(NA,dim=c(length(n),2,length(grid)),
         dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))

for (i in 1:length(n)){
  posterior.skpt<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
  }
  skpt.nc<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
  posterior.nc.skpt<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/skpt.nc
  }
  posterior.enth<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
  }
  enth.nc<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
  posterior.nc.enth<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/enth.nc
  }
  
  dens[i,"skpt",]<-posterior.nc.skpt(grid)*(-1)+spacing[i]
  dens[i,"enth",]<-posterior.nc.enth(grid)+spacing[i]
  loc[i,"skpt",]<-grid
  loc[i,"enth",]<-grid
  x[i,"skpt",]<-rep(as.character(i),length(grid))
  x[i,"enth",]<-rep(as.character(i),length(grid))
  m[i,"skpt",]<-rep("s",length(grid))
  m[i,"enth",]<-rep("e",length(grid))
}
#dens[1,"skpt",1]<-0 # 12/8/19
#dens[1,"enth",length(grid)]<-0 # 12/8/19
pdat=data.frame(dens=as.vector(dens),
                loc=as.vector(loc),
                x=as.factor(as.vector(x)),
                m=as.factor(as.vector(m)))

pdat2<-pdat[order(pdat$x,pdat$m),]

plot(0,type='n',ann=FALSE,
     xlim=c(-5,max(spacing)+5),
     ylim=c(0,1),
     xlab="",
     ylab="",
     axes=FALSE)

# plot(pdat2$dens[pdat2$x==1 & pdat2$m=="s"],
#      pdat2$loc[pdat2$x==1 & pdat2$m=="s"],
#      type='l',
#      xlim=c(-5,max(spacing)+5),
#      xlab="",
#      ylab="",
#      axes=FALSE)
box()
abline(h=seq(0,1,by=0.2),col='grey')

title(ylab="Posterior Probability of Response", line=3)
axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
axis(1,las=0,at=spacing,labels=rep("",length(spacing)))

for (i in seq(1,length(n))){
 for (j in c("e")){
#lines(pdat2$dens[pdat2$x==i & pdat2$m==j],
#     pdat2$loc[pdat2$x==i & pdat2$m==j])
   
   polygon(y=c(spacing[i],pdat2$loc[pdat2$x==i & pdat2$m==j],1,1),
           x=c(spacing[i],pdat2$dens[pdat2$x==i & pdat2$m==j],pdat2$dens[pdat2$x==i & pdat2$m==j][x.len],spacing[i]),
           col="lightgrey", border = NA)
 }
}
for (i in seq(1,length(n))){
  for (j in c("s")){
#    lines(pdat2$dens[pdat2$x==i & pdat2$m==j],
#          pdat2$loc[pdat2$x==i & pdat2$m==j])
    
    polygon(y=c(spacing[i],pdat2$loc[pdat2$x==i & pdat2$m==j],1,1),
            x=c(spacing[i],pdat2$dens[pdat2$x==i & pdat2$m==j],pdat2$dens[pdat2$x==i & pdat2$m==j][x.len],spacing[i]),
            col="darkgrey", border = NA)
  }
}

# Remove posterior means 7/21/2021
# for (i in 1:length(n)){
# points(spacing[i],lower[i],pch="-",cex = 2)
# points(spacing[i],upper[i],pch="-",cex = 2)
# points(spacing[i],mean[i],pch=20)
# lines(rep(spacing[i],2),c(lower[i],upper[i]),lwd=1)
# }
# 
# for (i in 2:length(n)){
# lines(c(spacing[i-1],spacing[i]),c(mean[i-1],mean[i]),lwd=1)
#   }

row<-1
mtext(text="Enrolled",side=1,line=row,at=label.x,adj=0)
for (i in 1:length(n)){
  mtext(text=n[i]+miss[i],side=1,line=row,at=spacing[i],adj=0)
}


row<-2
mtext(text="Completed",side=1,line=row,at=label.x,adj=0)
for (i in 1:length(n)){
  mtext(text=n[i],side=1,line=row,at=spacing[i],adj=0)
}

row<-3
mtext(text="Response Rate",side=1,line=row,at=label.x,adj=0)
for (i in 1){
  mtext(text="n/a",
        side=1,line=row,at=spacing[i],adj=0)
}
for (i in 2:length(n)){
  mtext(text=paste0(format(round(y1[i]/n[i]*100,digits=1),nsmall=1),"%"),
        side=1,line=row,at=spacing[i],adj=0)
}

row<-4
mtext(#text=as.expression(bquote(P(theta>(theta[0]+theta[1])/2 *"|"* bold(D) , pi[S]))),
      text="Efficacy Criteria",
      side=1,line=row,at=label.x,adj=0)
for (i in 1:length(n)){
  mtext(text=paste0(format(round(1 - efficacy[i],digits=3),nsmall=3)),side=1,line=row,at=spacing[i],adj=0)
}

row<-5
mtext(#text=as.expression(bquote(P(theta<theta[0] *"|"* bold(D) , pi[E]))),
      text="Futility Criteria",
      side=1,line=row,at=label.x,adj=0)
for (i in 1:length(n)){
  mtext(text=paste0(format(round(futility[i],digits=3),nsmall=3)),side=1,line=row,at=spacing[i],adj=0)
}

#lines(rep(spacing[2]/2,2),c(0,1),lty="dashed")
#lines(rep((spacing[5]+spacing[4])/2,2),c(0,1),lty="dashed")