setwd("/Users/evankwiatkowski/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")

dat <- read.csv("Table0_merged.csv")

head(dat)
table(dat$p.IP)

dat$success <- dat$eff.prob.initial >= 0.975

dat$pm.rd.f2 <- (dat$y1.IP.f / (dat$y1.IP.f + dat$y0.IP.f)) - (dat$y1.PC.f / (dat$y1.PC.f + dat$y0.PC.f))





final <- aggregate(dat, list(dat$p.IP), mean)
final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
final$mean.ss2 <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss

final$mean.ss
final$mean.ss2
final$pm.rd.f
final$pm.rd.f2

final$success
# subs <- dat[dat$p.IP == 0.54, !names(dat) %in% c("X.1", "X")]
# colMeans(subs)
# mean(subs$y1.IP.f + subs$y0.IP.f + subs$y1.PC.f + subs$y0.PC.f) # avg ss


# x.seq <- final$p.IP - final$p.PC
# x.seq
# plot(x.seq, x.seq, 
#      type = 'l', 
#      lty = 'dashed',
#      # xlim = c(-0.03, 0.15), 
#      # ylim = c(-0.03, 0.15),
#      xlab = "True Treatment Effect Value",
#      ylab = "Average Posterior Mean")
# lines(x.seq, final$pm.rd.f, lwd = 3)
# points(x.seq, final$pm.rd.f, pch = 8)
# lines(x.seq, final$pm.ni.f, col = "green", lwd = 3)
# lines(x.seq, final$pm.skpt.enth.ni.f, lwd = 3)
# points(x.seq, final$pm.skpt.enth.ni.f, pch = 16)
# 
# 
# 
# # lines(x.seq, final$pm.rd.f, type = 'l', lwd = 3)
# # lines(x.seq, final$pm.skpt.enth.ni.tilde.f)
# 
# # lines(x.seq, final$pm.skpt.enth.f)
# 
# lines(x.seq, final$pm.skpt.f, col = 'blue', lwd = 3)
# # abline(v = 0)
# lines(x.seq, final$pm.enth.f, col = 'red', lwd = 3)
# # abline(v = 0.12)
# # lines(x.seq, final$pm.50.50.f)
# # abline(v = 0.06)
# 
# legend("topleft",
#        c("MLE",
#          "Skpt",
#          "Enth",
#          "Non-informative",
#          # "50/50 mixture skpt/enth",
#          # "Adaptive skpt/enth",
#          "Adaptive skpt/enth/ni",
#          # "Adaptive skpt/enth/ni v2",
#          "Reference Line"),
#        lwd = 2,
#        lty = c(rep(1,5), 2),
#        col = c("black", "blue", "red", "green", "black", "black"),
#        pch = c(8, NA, NA, NA, 16, NA)
# )
# 
