setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/")

Table0 <- read.csv("../../output/Table0_merged.csv")
Table0 <- Table0[Table0$eff.mix.prob.x==10,]

par(mfrow=c(1,2))

x <- c(Table0$risk.diff.initial,Table0$risk.diff.final)
y <- c(Table0$box.skpt.initial,Table0$box.skpt.final)
plot(x, y, xlab = "Observed Response Difference", ylab = "Box's p-value", pch = 19, cex = 0.25)
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], lwd=2 )  

y <- c(Table0$box.enth.initial,Table0$box.enth.final)
points(x,y, col=2, pch = 19, cex = 0.25)
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], col=2, lwd=2 )  

y <- c(Table0$box.ni.initial,Table0$box.ni.final)
points(x,y, col="blue", pch = 19, cex = 0.25)
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], col="blue", lwd=2 )  

abline(v=0, col="grey")
abline(v=0.12, col="grey")

omega.skpt <- Table0$box.skpt.initial
omega.enth <- Table0$box.enth.initial
omega.ni   <- Table0$box.ni.initial
omega.ni   <- pmax(omega.ni-pmax(omega.skpt,omega.enth),0)
omega.sum  <- omega.skpt+omega.enth+omega.ni
omega.skpt <- omega.skpt/omega.sum
omega.enth <- omega.enth/omega.sum
omega.ni   <- omega.ni/omega.sum

x <- c(Table0$risk.diff.initial)
y <- omega.skpt
plot(x,y, pch = 19, cex = 0.25, ylim = c(0,1), xlab = "Observed Response Difference", ylab = "Mixture Weights")
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], lwd=2 )  

y <- omega.enth
points(x,y, col=2, pch = 19, cex = 0.25)
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], lwd=2, col = 2)

y <- omega.ni
points(x,y, col="blue", pch = 19, cex = 0.25)
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], lwd=2, col = "blue")  

## Used for 7/2/20 plot
# y1 <- c(Table0$box.skpt.initial,Table0$box.skpt.final)
# y2 <- c(Table0$box.enth.initial,Table0$box.enth.final)
# y  <- y2 - y1
# plot(x, y, xlab = "Observed Risk Difference", ylab = "Difference in Box's p-value")
# model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
# myPredict <- predict( model ) 
# ix <- sort(x,index.return=T)$ix
# lines(x[ix], myPredict[ix], lwd=2 )  
# abline(v=0, col="grey")
# abline(v=0.12, col="grey")
# abline(h=0)