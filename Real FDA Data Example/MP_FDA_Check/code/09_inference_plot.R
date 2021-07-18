rm(list	=	ls())
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/Real FDA Data Example/MP_FDA_Check/output/Table0")

# Merge output files and summarize
overall_list <- list.files()
for (i in 1:length(overall_list)){
  temp <- read.csv(file = paste0(overall_list[i]), header = TRUE, sep = ",")[,-1]
  if (i==1) {
    overall <- rbind(temp)
  } else {
    overall <- rbind(overall, temp)
  }
}

final <- overall
final$theta <- (final$y1.IP.f / (final$y1.IP.f + final$y0.IP.f)) - (final$y1.PC.f / (final$y1.PC.f + final$y0.PC.f))
plot(final$theta, final$pm.rd.f, type = 'l')
lines(final$theta, final$pm.skpt.f)
lines(final$theta, final$pm.enth.f)
lines(final$theta, final$pm.ni.f)
lines(final$theta, final$pm.50.50.f)
lines(final$theta, final$pm.skpt.enth.f)
lines(final$theta, final$pm.skpt.enth.ni.f)
lines(final$theta, final$pm.skpt.enth.ni.tilde.f)
