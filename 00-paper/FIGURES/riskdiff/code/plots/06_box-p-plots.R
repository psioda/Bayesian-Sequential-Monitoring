Table0 <- read.csv("../../output/Table0_merged.csv")
Table0 <- Table0[Table0$eff.mix.prob.x==10,]

par(mfrow=c(1,2))

x <- c(Table0$risk.diff.initial,Table0$risk.diff.final)
y <- c(Table0$box.skpt.initial,Table0$box.skpt.final)
plot(x, y, xlab = "Observed Response Difference", ylab = "Box's p-value")
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], lwd=2 )  

y <- c(Table0$box.enth.initial,Table0$box.enth.final)
points(x,y, col=2)
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], col=2, lwd=2 )  
abline(v=0, col="grey")
abline(v=0.12, col="grey")

y1 <- c(Table0$box.skpt.initial,Table0$box.skpt.final)
y2 <- c(Table0$box.enth.initial,Table0$box.enth.final)
y  <- y2 - y1


plot(x, y, xlab = "Observed Risk Difference", ylab = "Difference in Box's p-value")
model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
myPredict <- predict( model ) 
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix], lwd=2 )  
abline(v=0, col="grey")
abline(v=0.12, col="grey")
abline(h=0)
