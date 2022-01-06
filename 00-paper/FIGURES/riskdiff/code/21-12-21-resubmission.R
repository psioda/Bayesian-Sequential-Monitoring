rm(list = ls())
dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-4B.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-1-01-05-22.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-05-01-05-22.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-real-NI-01-04-22.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-01-04-22.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-01-04-22.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-T1E-12-28-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-NI-power-12-28-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-T1E-12-28-21.csv")
# dat <- read.csv("/Users/ek50/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output/Table0_merged-101-power-12-28-21.csv")

summary(dat$eff.prob.final)
summary(dat$min.ss)
summary(dat$max.ss)
# summary(dat$y1.IP.initial + dat$y0.IP.initial + dat$y1.PC.initial + dat$y0.PC.initial)
# summary(dat$y1.IP.final + dat$y0.IP.final + dat$y1.PC.final + dat$y0.PC.final)

# hist(dat[dat$p.IP == 0.39, "eff.prob.final"], freq = F, breaks = 10)
# abline(h = 1)

dat$reject <- (dat$eff.prob.final > 0.975)
dat$reject <- (dat$eff.prob.final > 0.887)
# dat$reject <- (dat$eff.prob.final > 0.787)
# dat$reject <- (dat$eff.prob.final > 0.72)

dat3 <- aggregate(dat$reject, list(dat$p.IP), FUN = mean) 
dat3$Group.1 <- dat3$Group.1 - 0.39
dat3

width.scale <- 9
output_png <- TRUE
if(output_png){
  png('7B.png',
      width = 450*width.scale,
      height = 300*width.scale,
      pointsize=16,
      res=300)
}
# par(mar=c(5.1 + 5, 4.1 + 2, 2.1, 2.1 + 0.5)) #c(bottom, left, top, right)

plot(dat1$Group.1, 
     dat1$x, 
     type = 'l', 
     lty = 1, 
     lwd = 1,
     xlim = c(-0.1, 0.3), 
     xlab = "Risk Difference", 
     ylab = "Null Hypothesis Rejection Rate", 
     xaxt = "n",
     yaxt = "n")
abline(h = seq(0, 1, by = 0.1),
       col = 'grey')
axis(1,
     at = seq(-.1, .3, by = .05),
     labels = seq(-.1, .3, by = .05))
axis(2, 
     at = seq(0, 1, by = 0.1), 
     labels = seq(0, 1, by = 0.1)) 
points(dat1[dat1$Group.1==0, "Group.1"], dat1[dat1$Group.1 == 0, "x"])
# abline(h = 0.025, col = 'gray')
lines(dat2$Group.1, dat2$x, lty = 2, lwd = 3)
points(dat2[dat2$Group.1==0, "Group.1"], dat2[dat2$Group.1 == 0, "x"])
lines(dat3$Group.1, dat3$x, lty = 3, lwd = 3)
legend(x="topleft", 
       title = 'Prior for Treatment Efficacy',
       legend = c("Adaptive Monitoring Prior", "Non-Informative", "Non-Informative with Modified Critical Value"),
       lty = c(1, 2, 3),
       lwd = 3
)
mtext("(B)", side = 2, line = 3, at = 1, las = 1)

if(output_png){dev.off()}




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################




library(fmsb)
res <- riskdifference(dat$y1.IP.final,
                      dat$y1.PC.final,
                      dat$y1.IP.final + dat$y0.IP.final,
                      dat$y1.PC.final + dat$y0.PC.final)
summary(res$p.value)
mean(res$p.value < 0.1)
mean(res$p.value < 0.05)
mean(res$p.value < 0.025)

dat$conf.lower <- res$conf.int[1:nrow(dat)]
dat$conf.upper <- res$conf.int[(nrow(dat) + 1):(2 * nrow(dat))]


# dat$freq.rej <- (res$p.value < 0.025)
dat$freq.rej <- (dat$conf.lower > 0) #| (dat$conf.upper < 0)
aggregate(dat$conf.lower, list(dat$p.IP), FUN = mean) 
aggregate(dat$freq.rej, list(dat$p.IP), FUN = mean) 



res <- riskdifference(dat$y1.PC.final[1:3],
                      dat$y1.IP.final[1:3],
                      dat$y1.PC.final[1:3] + dat$y0.PC.final[1:3],
                      dat$y1.IP.final[1:3] + dat$y0.IP.final[1:3])
res$p.value
res$conf.int





# merge1 <- aggregate(dat$reject, list(dat$p.IP), FUN = mean)
merge2 <- aggregate(dat$reject, list(dat$p.IP), FUN = mean)
merge1$y <- merge2$x
merge1$z <- merge1$y - merge1$x
merge1

hist(dat[dat$p.IP == 0.39, "eff.prob.final"], freq = F)
abline(h = 1)







colnames(merge2) <- c("g", "y")
merge <- cbind(merge1, merge2)

plot(merge$g, merge$x, type = 'l', pch = 10)
lines(merge$g, merge$y)


mean(dat$eff.prob.final > 0.9)
mean(dat$eff.prob.final > 0.95)
new.alpha <- 0.2640385
mean(dat$eff.prob.final >= 1 - new.alpha)
# mean(dat$eff.prob.final > 0.975)
# mean(dat$eff.prob.final > 0.99)

# new.alpha <- 0.377

# mean(dat$eff.prob.final >= 0.9796)
# mean(dat$eff.prob.final >= 0.819)




# mean(res$p.value < new.alpha)

power.prop.test(n = 266, p1 = .39, p2 = .51, sig.level = 0.025, alternative = "two.sided") 
power.prop.test(n = 269, p1 = .39, p2 = .51, sig.level = 0.025, alternative = "one.sided") 
power.prop.test(n = 210, p1 = .39, p2 = .51, sig.level = 0.025, alternative = "one.sided") 

power.prop.test(n = 161, p1 = .39, p2 = .51, sig.level = 0.05, alternative = "one.sided") 

## CHANGE IN PAPER TOO