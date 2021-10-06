rm(list = ls())
library(RColorBrewer)

setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")
dat <- read.csv("Table0_merged_2021-07-21-bayes-p-grid-90.csv")
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/2d heatmap bayes p")

head(dat)
summary(dat$eff.mix.prob.initial)
summary(dat$eff.mix.prob.final)
summary(dat$box.enth.initial)
summary(dat$box.skpt.initial)

output_png <- TRUE
par(mar=c(4.1, 4.1, 1.1, 1.1)) # c(bottom, left, top, right))
width.scale<-8
# if(output_png){png('2dbayesp.png',width = 300*1.25*width.scale, height = 300*width.scale,pointsize=16,res=300)}
if(output_png){png('2dpostp.png',width = 300*1.25*width.scale, height = 300*width.scale,pointsize=16,res=300)}
par(mar=c(4.1, 4.1, 1.1, 1.1)) # c(bottom, left, top, right))


grid.skpt <- dat
plot(NULL,
     xlim = c(0, 52),
     ylim = c(0, 38),
     xlab="",
     ylab="",
     col="white",
     xaxt='n',
     yaxt='n',
     xaxs = "i", 
     yaxs = "i")

axis(1,
     at = 0:10*(52/10),
     labels = format(round(seq(0,1,by=0.1),digits=1),nsmall=1))
axis(2,
     at = 0:10*(38/10),
     labels = format(round(seq(0,1,by=0.1),digits=1),nsmall=1))

grid.skpt$x <- grid.skpt$y1.IP.initial
grid.skpt$y <- grid.skpt$y1.PC.initial
# grid.skpt$z <- grid.skpt$box.enth.initial
grid.skpt$z <- grid.skpt$eff.prob.initial

title(ylab=as.expression(bquote("Control Sample Proportion")),line=3)
title(xlab=as.expression(bquote("Treatment Sample Proportion")),line=3)

# cuts<-c(0,1E-11,1E-12,1E-10,1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1,10,20,30)


cuts <- c(0, 0.025, 0.5, 0.975, 1)
cuts.names <- c("<0.025", "0.025 - 0.5", "0.5 - 0.975", ">0.975")
# 
# cuts <- c(0, 0.1, 0.25, 0.75, 0.9, 1)
# cuts.names <- c("<0.1", "0.1 - 0.25", "0.25 - 0.75", "0.75 - 0.9", ">0.9")

# colors <- gray.colors(length(cuts)-1, start = 0.9, end = 0)
colors <- brewer.pal(length(cuts)-1,"Blues")

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

points(27, 16, pch = 19, col = "white")
text(27, 16, labels = round(grid.skpt[grid.skpt$x == 27 & grid.skpt$y == 16, "z"],3), col = "white", pos = 1)
library(RColorBrewer)
abline(a = c(0, 0), b = 38/52, lwd = 3, lty = 'dashed')
legend(x="topleft", 
       # title = "Box's p-value",
       title = "Efficacy Posterior Probability",
       legend=c(cuts.names, "Observed Data", "Equivalent Proportions"),
       # legend=c(cuts.names, "Observed Data"),
       fill = c(colors,'black', NA),
       # col = c(rep(NA, length(colors)), 'black'),
       density = c(rep(NA, length(colors)), 0, 0),
       lty = c(rep(NA, length(colors)), NA, 'dashed'),
       pch = c(rep(NA, length(colors)), 21, NA),
       lwd = c(rep(NA, length(colors)), 1, 3),
       border = c(rep('black', length(colors)), NA, NA),
       pt.bg = c(rep(NA, length(colors)), 'white', NA),
       # cex = c(rep(NA, length(colors)), 2)
       )


# ## FILL IN LAST SEGMENT
# polygon(c(outer.x,outer.x[1]),
#         c(outer.y,outer.y[1]),
#         col=colors[i], 
#         border = NA)
# mtext("(B)", side = 2, line = 3, at = 38, las = 1)


if(output_png){dev.off()}
