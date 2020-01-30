#par(mfrow=c(1,2))

# print the overall labels
#source("plots/plots_priors.R")
#mtext('Density Value', side = 1, outer = TRUE, line = 3)
#mtext('Density Value', side = 2, outer = TRUE, line = .3)

### Prior plots  ########################################
# png(filename=paste0("normal_priors",prior_indicator_skpt[idx],prior_indicator_enth[idx],".png"))
# source("plots/plots_priors.R")
# dev.off()

# Step 5: Plot Skeptical and enthuastic priors separately
## 7/19/19 match slide 14 from FDA presentation, shade in gray scale, set 350 DPI

#par(ask=FALSE)
#par(oma = c(1.25, 1.25, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
#par(mar=c(5.1-1,4.1,4.1-4,2.1),mai = c(0.3, 0.3, 0.1, 0.1)) #bottom, left, top, and right.
#par(mfrow = c(1, 2)) 
x<-seq(0,1,by=0.01)


plot(x,prior.nc.skpt(x),type="l",
     xlab="",
     ylab="",
     main="",
     xaxt="n",
     yaxt="n",
     ylim=c(0,6)) # 20-01-02
title(ylab="Density Value", line=1, #cex.lab=1.2
)
title(xlab="Response Probability",line=2)
axis(1,at=c(p.enth,p.skpt),
     labels=c(as.expression(bquote(theta[1])),as.expression(bquote(theta[0]))))
#title(xlab="Density Value",line=2)
#axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
polygon(c(0,x[x<=p.enth],p.enth),c(0,prior.nc.skpt(x)[x<=p.enth],0),col="lightgrey")
polygon(c(p.enth,x[x>=p.enth],1),c(0,prior.nc.skpt(x)[x>=p.enth],0),col="black")
segments(x0=p.skpt,y0=0,y1=prior.nc.skpt(p.skpt))

legend('topright',
       legend= c(as.expression(bquote(mu == theta[0])),
                #as.expression(bquote(alpha == .(sigma0.skpt))),
                #as.expression(bquote(beta == .(lambda0.skpt))),
                as.expression(bquote(P(theta> theta[1])==.(tail.skpt)))),
       #inset=c(0.01,0.01)
       )
#text(x=0.7,y=5,labels=bquote(P(theta>.(p.enth))==.(tail.skpt)))
#text(x=0.7,y=4,labels=bquote(mu == .(p.skpt)),adj=0)
#text(x=0.7,y=3,labels=bquote(alpha == .(sigma0.skpt)),adj=0)
#text(x=0.7,y=2,labels=bquote(beta == .(lambda0.skpt)),adj=0)