x<-seq(0,1,by=0.01)

plot(x,prior.nc.enth(x),type="l",
     xlab="",
     ylab="",
     main="",
     #xaxt="n",
     #yaxt="n",
     ylim=c(0,3)) # 20-01-02
#axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
title(ylab="Response Probability", line=2, #cex.lab=1.2
)
title(xlab="Density Value",line=2)

polygon(c(x[x<=p.skpt],p.skpt),c(prior.nc.enth(x)[x<=p.skpt],0),col="black")
polygon(c(p.skpt,x[x>=p.skpt],1),c(0,prior.nc.enth(x)[x>=p.skpt],0),col="lightgrey")
segments(x0=p.enth,y0=0,y1=prior.nc.enth(p.enth))

legend('topleft',
       legend= c(as.expression(bquote(P(theta<.(p.skpt))==.(tail.enth))), 
                 as.expression(bquote(mu == .(p.enth))),
                 as.expression(bquote(alpha == .(sigma0.enth))),
                 as.expression(bquote(beta == .(lambda0.enth)))))

#text(x=0.1,y=2.5,labels=bquote(P(theta<.(p.skpt))==.(tail.enth)))
#text(x=0.1,y=2,labels=bquote(mu == .(p.enth)),adj=0)
#text(x=0.1,y=1.5,labels=bquote(alpha == .(sigma0.enth)),adj=0)
#text(x=0.1,y=1,labels=bquote(beta == .(lambda0.enth)),adj=0)
#text(x=0.1,y=1,labels=bquote(P(theta>.(lambda0.enth))==.(tail_skpt))))
