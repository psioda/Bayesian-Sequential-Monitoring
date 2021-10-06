width.scale<-5
par(oma=c(0, 0, 0, 0))
par(mar=c(0, 0, 0, 0))
png('enth_aug12.png',width = 300*2*width.scale, height = 300*width.scale,pointsize=16,res=300)
par(oma=c(0, 0, 0, 0))
par(mar=c(0, 0, 0, 0))
# assemble final prior
enth.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
  exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
    exp(-(abs(y - mu)/enth.alpha0)^enth.beta0)/(2*enth.alpha0*gamma(1/enth.beta0)/enth.beta0)/
    (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
       pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
    (pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
       pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
}
enth.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
  exp(-(abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0)/(2*enth.rd.alpha0*gamma(1/enth.rd.beta0)/enth.rd.beta0)*
    exp(-(abs(y - mu)/enth.alpha0)^enth.beta0)/(2*enth.alpha0*gamma(1/enth.beta0)/enth.beta0)/
    (pgnorm(q = 1,  mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0) -
       pgnorm(q = -1, mu = delta.enth, alpha = enth.rd.alpha0, beta = enth.rd.beta0))/
    (pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
       pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
}

x.len <- 101
y.len <- 101
x <- seq(-1, 1, length = x.len)
y <- seq(0,  1, length = y.len)
grid <- expand.grid(x = x, y = y)

grid.enth <- grid
for (i in 1:nrow(grid)){
  if (grid.enth$x[i] >= 0 && (grid.enth$x[i] < (1 - grid.enth$y[i]))) {
    grid.enth$z[i] <- enth.prior.1(grid.enth$x[i], grid.enth$y[i])
  } else if (grid.enth$x[i] <= 0 && (grid.enth$x[i] > -grid.enth$y[i])) {
    grid.enth$z[i] <- enth.prior.2(grid.enth$x[i], grid.enth$y[i])
  } else {
    grid.enth$z[i] <- 0
  }
}

par(mar=c(3.1, 3.1, 0.1, 0.1)) # c(bottom, left, top, right))

plot(grid.enth$x,grid.enth$y,
     xlab="",
     ylab="",
     col="white",
     xaxt='n',
     yaxt='n',
     xlim = c(-1,1),
     ylim = c(0,1))

title(ylab=as.expression(bquote("Response Probability")), line = 2)
title(xlab=as.expression(bquote("Risk Difference")), line = 2)

cuts<-c(0,1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)
# colors<-gray.colors(length(cuts)-1, start = 0.9, end = 0)
colors = rev(hcl.colors(length(cuts)-1, palette = "blues", alpha = NULL, rev = FALSE, fixup = TRUE))


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

axis(1,at=c(delta.skpt,delta.enth,-1,1),
     labels=c(as.expression(bquote(theta[0])),as.expression(bquote(theta[1])),-1,1))

axis(2,at=c(mu,0,1),
     labels=c(as.expression(bquote(eta[0])),0,1))

## FILL IN LAST SEGMENT
polygon(c(outer.x,outer.x[1]),
        c(outer.y,outer.y[1]),
        col=colors[i], 
        border = NA)

dev.off()

