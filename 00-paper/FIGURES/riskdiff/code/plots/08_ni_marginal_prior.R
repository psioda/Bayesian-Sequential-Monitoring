output_png <- FALSE

mu0.enth   <- delta.intr
tail.skpt  <- 0.025  # tail probabilities for priors
tail.enth  <- 0.025

# closed form
ni.rd.prior <- function(x){
  dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
    (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
       pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))
}

width.scale <- 6
if(output_png){png('figure1d.png', width = 300*width.scale,  height = 300*width.scale, pointsize=12, res=300)}
scale         <-  1.5

xmin  <-  -1
xmax  <-  1
ymax  <-  6
x <- seq(xmin, 
         xmax, 
         by=0.005)
plot(x, ni.rd.prior(x), type="l", 
     xlab="", 
     ylab="", 
     main="", 
     xaxt="n", 
     yaxt="n", 
     xlim=c(xmin, xmax), 
     ylim=c(0, ymax)) # 20-01-02
#axis(2, at=c(0, 1, 2, 3), labels=c(0, 1, 2, 3))
axis(1, at = seq(-.48,0.60,by=0.18), labels = round(seq(-.48,0.60,length=7),2))
# axis(1, at=c(delta.enth, delta.skpt, (delta.enth+delta.skpt)/2, (3*delta.skpt-delta.enth)/2), 
#      labels=c(as.expression(bquote(theta[1])), 
#               as.expression(bquote(theta[0])), 
#               as.expression(bquote((theta[0]+theta[1])/2)),
#               as.expression(bquote((3*theta[0]-theta[1])/2))))
title(ylab="Density Value",  line=1)
title(xlab="Response Probability", line=2)
#title(xlab="Density Value", line=2)

# polygon(c(x[x<=mu0.enth-(delta.enth-delta.skpt)], mu0.enth-(delta.enth-delta.skpt)), 
#         c(ni.rd.prior(x)[x<=mu0.enth-(delta.enth-delta.skpt)], 0), col="black")
# 
# polygon(c(mu0.enth-(delta.enth-delta.skpt), x[x>=mu0.enth-(delta.enth-delta.skpt) & x<=delta.skpt], delta.skpt), 
#         c(0, ni.rd.prior(x)[x>=mu0.enth-(delta.enth-delta.skpt) & x<=delta.skpt], 0), col="lightgrey")

segments(x0=delta.ni.enth, y0=0, y1=ni.rd.prior(delta.ni.enth))
segments(x0=delta.ni.skpt, y0=0, y1=ni.rd.prior(delta.ni.skpt))
segments(x0=delta.ni.intr, y0=0, y1=ni.rd.prior(delta.ni.intr))

# legend("top", 
#        legend= c(
#          as.expression(bquote(mode(theta) == (theta[0]+theta[1])/2)), 
#          as.expression(bquote(P(theta< (3*theta[0]-theta[1])/2)==.(tail.enth))), 
#          as.expression(bquote(P(theta %in% ((3*theta[0]-theta[1])/2*", "*theta[0])==.(round((pnorm(qnorm(tail.enth)/2)-tail.enth)*scale, 3)))))#, 
#          #as.expression(bquote(GN(mu==theta[1], alpha==.(sigma0.enth), beta==.(lambda0.enth))))
#        ))

mtext("(D)", side=2, line=1, at=ymax, las=1)
if(output_png){dev.off()}