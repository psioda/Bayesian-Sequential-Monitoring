# layout(matrix(c(2,2,4,
#                 1,1,3,
#                 1,1,3),
#               nrow=3,
#               byrow = TRUE),
#        widths=c(2,2,1),
#        heights=c(1,2,2))
# #layout.show(n=4)

par(mar=c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))
#par(mar=c(2.1, 2.1, 0, 0)) # c(bottom, left, top, right))

plot(grid.enth$x,grid.enth$y,xlab="",ylab="",col="white",xaxt='n',yaxt='n')
#axis(1,at=mu,labels=c(as.expression(bquote(mu))))
#axis(2,at=mu+delta.enth,labels=c(as.expression(bquote(mu + delta[E]))))

title(ylab=as.expression(bquote("Response Probability"~theta[1])),line=1)
title(xlab=as.expression(bquote("Response Probability"~theta[0])),line=1)

cuts<-c(1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)
colors<-gray.colors(length(cuts)-1, start = 1, end = 0)

for (i in 1:length(cuts)-1){
  outer.xy <- grid.enth[grid.enth$z>cuts[i+1],c("x","y")]
  inner.xy <- grid.enth[grid.enth$z>cuts[i],c("x","y")]
  outer.x<-outer.xy[chull(outer.xy),"x"]
  outer.y<-outer.xy[chull(outer.xy),"y"]
  inner.x<-inner.xy[chull(inner.xy),"x"]
  inner.y<-inner.xy[chull(inner.xy),"y"]
  polygon(c(outer.x,outer.x[1],inner.x,inner.x[1]),
          c(outer.y,outer.y[1],inner.y,inner.y[1]),
          col=colors[i], 
          border = NA)
}
## FILL IN LAST SEGMENT
polygon(c(outer.x,outer.x[1]),
        c(outer.y,outer.y[1]),
        col=colors[i], 
        border = NA)

## Add lines
lines(c(0,1-delta.enth),c(delta.enth,1))
lines(c(0,1),c(0,1))
#lines(c(0,mu),c(mu+delta.enth,mu+delta.enth)) # horizontal line
#lines(c(mu,mu),c(0,mu+delta.enth)) # vertical line
#par(mar=c(5.1, 4.1, 4.1, 2.1)/2)
#    oma=c(0,0,0,0)) # c(bottom, left, top, right))

legend('bottomright',
       legend= c(#as.expression(bquote(mu == .(p.enth))),
         #as.expression(bquote(delta[S] == .(delta.skpt)~","~delta[E] == .(delta.enth))),
         #as.expression(bquote(mu[0] == .(mu)*","~alpha[0] == .(sigma0.placebo)*","~beta[0] == .(lambda0.placebo))),
         #as.expression(bquote(delta == delta[E]*","~alpha[1] == .(sigma0.enth)*","~beta[1] == .(lambda0.enth))),
         as.expression(bquote(P(theta[1]<theta[0]+delta[S])==.(1-sig.eff)))))