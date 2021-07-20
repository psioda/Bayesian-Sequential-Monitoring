##################################
# Enrollment distribution
# Evan Kwiatkowski, Feb 2020
##################################

if (is.na(p.IP)){
  dat <- read.csv("../data/datasummary.csv")
  
  group             <- dat$group
  group[group == 0] <- "PC"
  group[group == 1] <- "IP"
  enr.times.all     <- dat$enr.times.all
  outcome.times.all <- dat$outcome.times.all
  responses         <- dat$responses
  
  enr.times.PC     <- enr.times.all[group == "PC"]
  outcome.times.PC <- outcome.times.all[group == "PC"]
  responses.PC     <- responses[group == "PC"]
  
  enr.times.IP     <- enr.times.all[group == "IP"]
  outcome.times.IP <- outcome.times.all[group == "IP"]
  responses.IP     <- responses[group == "IP"]
} else {
  group <- c(sample(c(rep("PC", 4),
                      rep("IP", 20))),
             sample(c(rep("PC", 38),
                      rep("IP", 38))))
  enr.times.all     <- seq(1:100) * 17.2
  outcome.times.all <- enr.times.all + 365.25
  
  enr.times.PC     <- enr.times.all[group == "PC"]
  outcome.times.PC <- outcome.times.all[group == "PC"]
  responses.PC     <- rbinom(n = 42, size = 1, prob = p.PC)
  
  enr.times.IP     <- enr.times.all[group == "IP"]
  outcome.times.IP <- outcome.times.all[group == "IP"]
  responses.IP     <- rbinom(n = 58, size = 1, prob = p.IP)
}