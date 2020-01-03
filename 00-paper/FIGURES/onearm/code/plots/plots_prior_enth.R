x<-seq(0,1,by=0.01)

plot(x,prior.nc.enth(x),type="l",
     xlab="",
     ylab="",
     main="",
     #xaxt="n",
     #yaxt="n",
     ylim=c(0,3)) # 20-01-02
#axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
polygon(c(x[x<=p.skpt],p.skpt),c(prior.nc.enth(x)[x<=p.skpt],0),col="black")
polygon(c(p.skpt,x[x>=p.skpt],1),c(0,prior.nc.enth(x)[x>=p.skpt],0),col="lightgrey")
segments(x0=p.enth,y0=0,y1=prior.nc.enth(p.enth))