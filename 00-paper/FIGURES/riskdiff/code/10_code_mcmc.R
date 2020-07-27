rm(list=ls())

library(LearnBayes)
library(MCMCpack)
library(coda)
library(gnorm)
library(foreach)
library(doParallel)
registerDoParallel(detectCores())
getDoParWorkers()

load(file = 'args_model.RData')

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# source("10A_code_mcmc_preprocessing.R")
Table0_merged                <- read.csv("../output/Table0_merged_071620.csv")
Table0_merged                <- Table0_merged[Table0_merged$eff.mix.prob==10,]

names.merged                 <- c("y1.PC","y0.PC","y1.IP","y0.IP", "box.skpt", "box.enth")
names.initial                <- c("y1.PC.initial","y0.PC.initial","y1.IP.initial","y0.IP.initial", "box.skpt.initial","box.enth.initial")
Table0_merged.initial        <- Table0_merged[!duplicated(Table0_merged[names.initial]), names.initial]
names(Table0_merged.initial) <- names.merged

names.final                  <- c("y1.PC.final","y0.PC.final","y1.IP.final","y0.IP.final", "box.skpt.final","box.enth.final")
Table0_merged.final          <- Table0_merged[!duplicated(Table0_merged[names.final]), names.final]
names(Table0_merged.final)   <- names.merged

Table0_merged                <- rbind(Table0_merged.initial, Table0_merged.final)
Table0_merged                <- Table0_merged[!duplicated(Table0_merged),]

Table0_ni_unique             <- read.csv("../output/Table0_ni_unique_072220.csv")

Table0_final                 <- merge(Table0_ni_unique, Table0_merged, by = c("y1.PC","y0.PC","y1.IP","y0.IP"))
Table0_final$box.ni          <- Table0_final$box.ni.final
Table0_final                 <- Table0_final[c("y1.PC","y0.PC","y1.IP","y0.IP", "box.ni", "box.skpt", "box.enth")]
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


## mcmc function
logpost <- function(theta, data) {
  x          <- theta[1] # rd
  y          <- theta[2] # rd
  y11        <- data[["y1.IP"]]
  y10        <- data[["y0.IP"]]
  y01        <- data[["y1.PC"]]
  y00        <- data[["y0.PC"]]
  
  omega.skpt <- data[["box.skpt"]]
  omega.enth <- data[["box.enth"]]
  omega.ni   <- data[["box.ni"]]
  omega.ni   <- pmax(omega.ni-pmax(omega.skpt,omega.enth),0)
  omega.sum  <- omega.skpt + omega.enth + omega.ni
  omega.skpt <- omega.skpt/omega.sum
  omega.enth <- omega.enth/omega.sum
  omega.ni   <- omega.ni/omega.sum
  
  llike    <-  dbinom(y11, y11 + y10, x + y, log = TRUE) +
               dbinom(y01, y01 + y00, y,     log = TRUE)
  
  prior    <-  log(
    omega.skpt*
      (dgnorm(x, delta.skpt,    skpt.rd.alpha0, skpt.rd.beta0)/(pgnorm(q = 1,    delta.skpt,    skpt.rd.alpha0, skpt.rd.beta0) - pgnorm(q = -1, delta.skpt,    skpt.rd.alpha0, skpt.rd.beta0))*
         ((x>=0)*dgnorm(y, mu,            skpt.alpha0,    skpt.beta0)   /(pgnorm(q = 1-x,  mu,            skpt.alpha0,    skpt.beta0) -    pgnorm(q = 0,  mu,            skpt.alpha0,    skpt.beta0)) +
            (x<0)*dgnorm(y, mu,            skpt.alpha0,    skpt.beta0)   /(pgnorm(q = 1,    mu,            skpt.alpha0,    skpt.beta0) -    pgnorm(q = -x, mu,            skpt.alpha0,    skpt.beta0)))) +
      omega.enth*
      (dgnorm(x, delta.enth,    enth.rd.alpha0, enth.rd.beta0)/(pgnorm(q = 1,    delta.enth,    enth.rd.alpha0, enth.rd.beta0) - pgnorm(q = -1, delta.enth,    enth.rd.alpha0, enth.rd.beta0))*
         ((x>=0)*dgnorm(y, mu,            enth.alpha0,    enth.beta0)   /(pgnorm(q = 1-x,  mu,            enth.alpha0,    enth.beta0) -    pgnorm(q = 0,  mu,            enth.alpha0,    enth.beta0)) +
            (x<0)*dgnorm(y, mu,            enth.alpha0,    enth.beta0)   /(pgnorm(q = 1,    mu,            enth.alpha0,    enth.beta0) -    pgnorm(q = -x, mu,            enth.alpha0,    enth.beta0)))) +
      omega.ni*
      (dgnorm(x, delta.ni.skpt, ni.rd.alpha0,   ni.rd.beta0)  /(pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0,   ni.rd.beta0) -   pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0,    ni.rd.beta0))*
         ((x>=0)*dgnorm(y, mu,            ni.alpha0,      ni.beta0)     /(pgnorm(q = 1-x,  mu,            ni.alpha0,      ni.beta0) -      pgnorm(q = 0,  mu,            ni.alpha0,       ni.beta0)) +
            (x<0)*dgnorm(y, mu,            ni.alpha0,      ni.beta0)     /(pgnorm(q = 1,    mu,            ni.alpha0,      ni.beta0) -      pgnorm(q = -x, mu,            ni.alpha0,       ni.beta0))))
  )
  val_post <- llike + prior
  return(val_post)
}

## data
burn_in          <- 10000  
nmcmc            <- 15000 + burn_in 
start            <- c(0, 0.5)

## analysis
start_time <- Sys.time()

x <- foreach (i = 1:nrow(Table0_final), .combine='c') %dopar% {
  # prior_dat_conflict(Table0$y1.IP.final[i],
  #                    Table0$y0.IP.final[i],
  #                    Table0$y1.PC.final[i],
  #                    Table0$y0.PC.final[i])
  data.all <- Table0_final[i,]
  laplace  = laplace(logpost, start, data.all) # change to data.all
  proposal = list(var = laplace$var, scale = 1)
  s        = rwmetrop(logpost, proposal, start, nmcmc, data.all)
  mcmc     = mcmc(s$par[-c(1:burn_in),])
  return(c(apply(mcmc, 2, mean), HPDinterval(mcmc)))
}

end_time  <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
cat("Started  ", as.character(start_time), "\n",
    "Finished ", as.character(end_time), "\n",
    "Time difference of ", diff_time, " ", attr(diff_time, "units"), "\n",
    sep = "")

## reformat final output
new.x                  <- matrix(x, ncol = 6,  byrow = TRUE)
colnames(new.x)        <- c("pm_rd","pm_pi0","lower_rd","upper_rd","lower_pi0","upper_pi0")
Table0_ni_unique_pm_cp <- cbind(Table0_final, new.x)

write.csv(Table0_ni_unique_pm_cp, "../output/Table0_ni_unique_pm_cp_072320.csv")