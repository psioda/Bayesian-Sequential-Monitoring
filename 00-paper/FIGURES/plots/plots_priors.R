# Step 5: Plot Skeptical and enthuastic priors separately
## 7/19/19 match slide 14 from FDA presentation, shade in gray scale, set 350 DPI

#par(ask=FALSE)
#par(oma = c(1.25, 1.25, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
#par(mar=c(5.1-1,4.1,4.1-4,2.1),mai = c(0.3, 0.3, 0.1, 0.1)) #bottom, left, top, and right.
#par(mfrow = c(1, 2)) 
x<-seq(0,1,by=0.01)


plot(x,prior.nc.skpt(x),type="l",
     xlab="",ylab="Density Value",
     main="",
     #xaxt="n",
     #yaxt="n",
     ylim=c(0,max(prior.nc.skpt(x))*1.1))
#axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
polygon(c(0,x[x<=p.enth],p.enth),c(0,prior.nc.skpt(x)[x<=p.enth],0),col="lightgrey")
polygon(c(p.enth,x[x>=p.enth],1),c(0,prior.nc.skpt(x)[x>=p.enth],0),col="black")
segments(x0=p.skpt,y0=0,y1=prior.nc.skpt(p.skpt))


plot(x,prior.nc.enth(x),type="l",
     xlab="Response Probability",
     ylab="",
     main="",
     #xaxt="n",
     #yaxt="n",
     ylim=c(0,max(prior.nc.enth(x))*1.1))
#axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
polygon(c(x[x<=p.skpt],p.skpt),c(prior.nc.enth(x)[x<=p.skpt],0),col="black")
polygon(c(p.skpt,x[x>=p.skpt],1),c(0,prior.nc.enth(x)[x>=p.skpt],0),col="lightgrey")
segments(x0=p.enth,y0=0,y1=prior.nc.enth(p.enth))