valentine <- function(index){
  
y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)

source("code_posteriors.R", local = TRUE)

print(paste0("PC.mle ", round(PC.mle, digits = 2),
             ", IP.mle ", round(IP.mle, digits = 2)))

print(paste0("skpt.lik ", round(skpt.lik, digits = 2),
             ", enth.lik ", round(enth.lik, digits = 2)))

print(paste0("eff.mix.prob ", round(eff.mix.prob, digits = 2)))

print(paste0("eff.mon.final ", round(inner[i, "eff.mon.final"], digits = 2),
             ", fut.mon.final ", round(inner[i, "fut.mon.final"], digits = 2)))

print("")
}