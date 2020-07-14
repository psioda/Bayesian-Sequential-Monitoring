# assemble final prior
ni.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
    (pgnorm(q = 1,     delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
       pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
    dgnorm(y,          mu,         ni.alpha0,    ni.beta0)/
    (pgnorm(q = 1 - x, mu,         ni.alpha0,    ni.beta0) -
       pgnorm(q = 0,   mu,         ni.alpha0,    ni.beta0))
}
ni.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
    (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
       pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
    dgnorm(y,         mu,         ni.alpha0,    ni.beta0)/
    (pgnorm(q = 1,    mu,         ni.alpha0,    ni.beta0) -
       pgnorm(q = -x, mu,         ni.alpha0,    ni.beta0))
}


x.len <- 101
y.len <- 101
x <- seq(-1, 1, length = x.len)
y <- seq(0,  1, length = y.len)
grid <- expand.grid(x = x, y = y)

grid.ni <- grid
for (i in 1:nrow(grid)){
  if (grid.ni$x[i] >= 0 && (grid.ni$x[i] < (1 - grid.ni$y[i]))) {
    grid.ni$z[i] <- ni.prior.1(grid.ni$x[i], grid.ni$y[i])
  } else if (grid.ni$x[i] <= 0 && (grid.ni$x[i] > -grid.ni$y[i])) {
    grid.ni$z[i] <- ni.prior.2(grid.ni$x[i], grid.ni$y[i])
  } else {
    grid.ni$z[i] <- 0
  }
}

par(mar=c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))

plot(grid.ni$x,grid.ni$y,
     xlab="",
     ylab="",
     col="white",
     xaxt='n',
     yaxt='n')

title(ylab=as.expression(bquote("Response Probability")), line = 3)
title(xlab=as.expression(bquote("Risk Difference")), line = 3)

cuts<-c(0,1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)
colors<-gray.colors(length(cuts)-1, start = 0.9, end = 0)

for (i in 1:length(cuts)-1){
  outer.xy <- grid.ni[grid.ni$z>cuts[i+1],c("x","y")]
  inner.xy <- grid.ni[grid.ni$z>cuts[i],c("x","y")]
  outer.x<-outer.xy[chull(outer.xy),"x"]
  outer.y<-outer.xy[chull(outer.xy),"y"]
  inner.x<-inner.xy[chull(inner.xy),"x"]
  inner.y<-inner.xy[chull(inner.xy),"y"]
  polygon(c(outer.x,outer.x[1],inner.x,inner.x[1]),
          c(outer.y,outer.y[1],inner.y,inner.y[1]),
          col=colors[i], 
          border = NA)
}

axis(1,at=c(delta.ni.skpt,delta.ni.enth,-1,1),
     labels=c(as.expression(bquote(theta[0])),as.expression(bquote(theta[1])),-1,1))

axis(2,at=c(mu,0,1),
     labels=c(as.expression(bquote(mu[0])),0,1))

## FILL IN LAST SEGMENT
polygon(c(outer.x,outer.x[1]),
        c(outer.y,outer.y[1]),
        col=colors[i], 
        border = NA)

