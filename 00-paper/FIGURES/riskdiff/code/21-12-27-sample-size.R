rm(list = ls())
library(fmsb)

SS <- 262
y1.PC.final <- rbinom(1E4, SS, prob = 0.39)
y1.IP.final <- rbinom(1E4, SS, prob = 0.39)

res <- riskdifference(y1.PC.final,
                      y1.IP.final,
                      SS,
                      SS)

summary(res$p.value)
# mean(res$p.value < 0.1)
# mean(res$p.value < 0.05)
mean(res$p.value < 0.025)

# https://clincalc.com/stats/samplesize.aspx