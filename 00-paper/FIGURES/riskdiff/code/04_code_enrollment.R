##################################
# Enrollment distribution
# Evan Kwiatkowski, Feb 2020
##################################
group <- c(sample(c(rep("PC", 4),
                    rep("IP", 20))),
           sample(c(rep("PC", 38),
                    rep("IP", 38))))
enr.times.all     <- seq(1:100) * 17.2
outcome.times.all <- enr.times.all + 52

enr.times.PC     <- enr.times.all[group == "PC"]
outcome.times.PC <- outcome.times.all[group == "PC"]
responses.PC     <- rbinom(n = 42, size = 1, prob = p.PC)

enr.times.IP     <- enr.times.all[group == "IP"]
outcome.times.IP <- outcome.times.all[group == "IP"]
responses.IP     <- rbinom(n= 58, size = 1, prob = p.IP)