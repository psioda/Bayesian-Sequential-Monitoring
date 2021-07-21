rm(list = ls())
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")
dat <- read.csv("Table0_merged_2021-07-21-bayes-p-grid.csv")
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/2d heatmap bayes p")

head(dat)
summary(dat$eff.mix.prob.initial)
summary(dat$eff.mix.prob.final)
summary(dat$box.enth.initial)
summary(dat$box.skpt.initial)

output_png <- FALSE
par(mar=c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))
width.scale<-6
if(output_png){png('2dbayesp.png',width = 300*2*width.scale, height = 300*width.scale,pointsize=16,res=300)}


grid.skpt <- dat
plot(NULL,
     xlim = c(0, 58),
     ylim = c(0, 42),
     xlab="",
     ylab="",
     col="white",
     xaxt='n',
     yaxt='n')

axis(1,
     at = c(0:58),
     labels = format(round((0:58)/58,digits=2),nsmall=2))
axis(2,
     at = c(0:42),
     labels = format(round((0:42)/42,digits=2),nsmall=2))

grid.skpt$x <- grid.skpt$y1.IP.initial
grid.skpt$y <- grid.skpt$y1.PC.initial
# grid.skpt$z <- grid.skpt$box.enth.initial
grid.skpt$z <- grid.skpt$eff.prob.initial

title(ylab=as.expression(bquote("Control Response Rate")),line=2)
title(xlab=as.expression(bquote("Treatment Response Rate")),line=3)

# cuts<-c(0,1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)

# cuts <- c(0,0.1,0.25,0.5,0.75,0.9,1)
cuts <- c(0, 0.1, 0.5, 0.975, 1)
colors <- gray.colors(length(cuts)-1, start = 0.9, end = 0)

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

points(31, 18, col = "white")

library(RColorBrewer)
legend(x="top", 
       legend=c("min", "ave", "max"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))

## FILL IN LAST SEGMENT
polygon(c(outer.x,outer.x[1]),
        c(outer.y,outer.y[1]),
        col=colors[i], 
        border = NA)



if(output_png){dev.off()}