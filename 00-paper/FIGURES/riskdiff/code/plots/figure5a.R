par(mar=c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))

plot(grid.skpt$x,grid.skpt$y,
     xlab="",
     ylab="",
     col="white",
     xaxt='n',
     yaxt='n')

title(ylab=as.expression(bquote("Response Probability"~theta[1])),line=1)
title(xlab=as.expression(bquote("Response Probability" ~theta[0])),line=1)

cuts<-c(1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)
colors<-gray.colors(length(cuts)-1, start = 1, end = 0)

for (i in 1:length(cuts)-1){
  outer.xy <- grid.skpt[grid.skpt$z>cuts[i+1],c("x","y")]
  inner.xy <- grid.skpt[grid.skpt$z>cuts[i],c("x","y")]
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

legend('topleft',
       legend= c(#as.expression(bquote(mu == .(p.enth))),
         #as.expression(bquote(delta[S] == .(delta.skpt)~","~delta[E] == .(delta.enth))),
         #as.expression(bquote(mu[0] == .(mu)*","~alpha[0] == .(alpha0.placebo)*","~beta[0] == .(beta0.placebo))),
         #as.expression(bquote(delta == delta[S]*","~alpha[1] == .(alpha0.skpt)*","~beta[1] == .(beta0.skpt))),
         as.expression(bquote(P(theta[1]>theta[0]+delta[E])==.(1-sig.eff)))))

