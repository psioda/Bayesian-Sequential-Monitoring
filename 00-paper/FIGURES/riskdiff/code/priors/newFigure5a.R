# assemble final prior
skpt.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
  exp(-(abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0)/(2*skpt.rd.alpha0*gamma(1/skpt.rd.beta0)/skpt.rd.beta0)*
    exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))/
    (pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
  exp(-(abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0)/(2*skpt.rd.alpha0*gamma(1/skpt.rd.beta0)/skpt.rd.beta0)*
    exp(-(abs(y - mu)/skpt.alpha0)^skpt.beta0)/(2*skpt.alpha0*gamma(1/skpt.beta0)/skpt.beta0)/
    (pgnorm(q = 1,  mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0) -
       pgnorm(q = -1, mu = delta.skpt, alpha = skpt.rd.alpha0, beta = skpt.rd.beta0))/
    (pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
       pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}

x.len <- 101
y.len <- 101
x <- seq(-1, 1, length = x.len)
y <- seq(0,  1, length = y.len)
grid <- expand.grid(x = x, y = y)

grid.skpt <- grid
for (i in 1:nrow(grid)){
  if (grid.skpt$x[i] >= 0 & (grid.skpt$y[i] < 1 - grid.skpt$x[i])) {
    grid.skpt$z[i] <- skpt.prior.1(grid.skpt$x[i], grid.skpt$y[i])
  } else if (grid.skpt$x[i] < 0 & (grid.skpt$y[i] > -grid.skpt$x[i])) {
    grid.skpt$z[i] <- skpt.prior.2(grid.skpt$x[i], grid.skpt$y[i])
  } else {
    grid.skpt$z[i] <- 0
  }
}

par(mar=c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))

plot(grid.skpt$x,grid.skpt$y,
     xlab="",
     ylab="",
     col="white",
     xaxt='n',
     yaxt='n')

title(ylab=as.expression(bquote("Response Probability")),line=3)
title(xlab=as.expression(bquote("Risk Difference")),line=3)

cuts<-c(0,1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)
colors<-gray.colors(length(cuts)-1, start = 0.9, end = 0)

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

#abline(h=c(0))
#abline(v=0)
## FILL IN LAST SEGMENT
polygon(c(outer.x,outer.x[1]),
        c(outer.y,outer.y[1]),
        col=colors[i], 
        border = NA)


axis(1,at=c(delta.skpt,delta.enth),
     labels=c(as.expression(bquote(theta[0])),as.expression(bquote(theta[1]))))

axis(2,at=c(mu),
     labels=c(as.expression(bquote(mu[0]))))