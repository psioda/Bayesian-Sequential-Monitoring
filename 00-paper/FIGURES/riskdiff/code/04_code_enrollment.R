##################################
# Enrollment distribution
# Evan Kwiatkowski, Feb 2020
##################################

if (is.na(p.IP)){
  
  ## BEGIN FAKE DATA FOR GRID OUTCOMES ##################
  # remember idx variable is available
  y0.seq <- 0:42
  y1.seq <- 0:58
  mat <- expand.grid(42, y0.seq, 58, y1.seq, 100, 100)
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