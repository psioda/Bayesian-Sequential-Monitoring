##################################
# Enrollment distribution
# Evan Kwiatkowski, Feb 2020
##################################

if (is.na(p.IP)){
  
  ## BEGIN FAKE DATA FOR GRID OUTCOMES ##################
  # remember idx variable is available
  y0.seq <- 0:38
  y1.seq <- 0:52
  mat <- expand.grid(38, y0.seq, 52, y1.seq, 92, 92)
  dat <- mat[idx, ]
  names(dat) <- c("nObs0", "yObs0", "nObs1", "yObs1", "targOutNum", "nFin")
  ## END FAKE DATA FOR GRID OUTCOMES ####################
  
  
  ## BEGIN REAL FDA FILE ################################
  # dat <- read.csv("../data/datasummary.csv")
  # 
  # group             <- dat$group
  # group[group == 0] <- "PC"
  # group[group == 1] <- "IP"
  # enr.times.all     <- dat$enr.times.all
  # outcome.times.all <- dat$outcome.times.all
  # responses         <- dat$responses
  # 
  # enr.times.PC     <- enr.times.all[group == "PC"]
  # outcome.times.PC <- outcome.times.all[group == "PC"]
  # responses.PC     <- responses[group == "PC"]
  # 
  # enr.times.IP     <- enr.times.all[group == "IP"]
  # outcome.times.IP <- outcome.times.all[group == "IP"]
  # responses.IP     <- responses[group == "IP"]
  ## END REAL FDA FILE ##################################
} else if (min.ss == 322){
  group <- sample(c(rep("PC", 161),
                    rep("IP", 161)))
  enr.times.all     <- seq(1:322) * 17.2
  outcome.times.all <- enr.times.all + 365.25
  
  enr.times.PC     <- enr.times.all[group == "PC"]
  outcome.times.PC <- outcome.times.all[group == "PC"]
  responses.PC     <- rbinom(n = 161, size = 1, prob = p.PC)
  
  enr.times.IP     <- enr.times.all[group == "IP"]
  outcome.times.IP <- outcome.times.all[group == "IP"]
  responses.IP     <- rbinom(n = 161, size = 1, prob = p.IP)

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