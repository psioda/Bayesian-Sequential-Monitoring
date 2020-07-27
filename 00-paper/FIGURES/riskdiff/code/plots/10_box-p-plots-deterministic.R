rm(list=ls())
library(pracma)
library(gnorm)
library(foreach)
library(doParallel)

registerDoParallel(detectCores())
getDoParWorkers()

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/")
load(file = '../args_model.RData') # loads all model information include prior parameters AND SETS SEED

# create data matrix
y1.IP        <- seq(0, 58)
y0.IP        <- 58 - y1.IP
Table0       <- data.frame(y1.IP, y0.IP)
Table0$y1.PC <- 16
Table0$y0.PC <- 26

prior_dat_conflict <- function(y1.IP, y0.IP, y1.PC, y0.PC){
  
  skpt.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  enth.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  ni.post.nc         <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)

  for (i in 1:(y1.IP + y0.IP + 1)){
    for (j in 1:(y1.PC + y0.PC + 1)){
      skpt.post.1   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
                  pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
            dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
            log(pgnorm(q = 1-x,  mu, skpt.alpha0, skpt.beta0) - 
                  pgnorm(q = 0,  mu, skpt.alpha0, skpt.beta0))
        )
      }
      skpt.post.2   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
                  pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
            dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    mu, skpt.alpha0, skpt.beta0) - 
                  pgnorm(q = -x, mu, skpt.alpha0, skpt.beta0))
        )
      }
      skpt.post.nc[i, j] <- integrate_debug(skpt.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
        integrate_debug(skpt.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
      
      enth.post.1   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
                  pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
            dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
            log(pgnorm(q = 1-x,  mu, enth.alpha0, enth.beta0) - 
                  pgnorm(q = 0,  mu, enth.alpha0, enth.beta0))
        )
      }
      enth.post.2   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
                  pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
            dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    mu, enth.alpha0, enth.beta0) - 
                  pgnorm(q = -x, mu, enth.alpha0, enth.beta0))
        )
      }
      enth.post.nc[i, j] <- integrate_debug(enth.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
        integrate_debug(enth.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
      
      ni.post.1 <- function(x, y){ # for x > 0 (theta > 0)
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
          log(pgnorm(q = 1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
             pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
          dgnorm(y,          mu,            ni.alpha0,    ni.beta0, log = TRUE) -
          log(pgnorm(q = 1-x,mu,            ni.alpha0,    ni.beta0) -
             pgnorm(q = 0,   mu,            ni.alpha0,    ni.beta0))
        )
      }
      ni.post.2 <- function(x, y){ # for x < 0 (theta < 0)
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
          log(pgnorm(q = 1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
             pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
          dgnorm(y,         mu,            ni.alpha0,    ni.beta0, log = TRUE) -
          log(pgnorm(q = 1, mu,            ni.alpha0,    ni.beta0) -
             pgnorm(q = -x, mu,            ni.alpha0,    ni.beta0))
        )
      }
      ni.post.nc[i, j] <- integrate_debug(ni.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                          integrate_debug(ni.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
    }
  }
  # prior data conflict for skeptical prior
  box.skpt <- sum(skpt.post.nc[skpt.post.nc <= skpt.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(skpt.post.nc)))
  print(paste0("Skeptical prior compatibility: ", box.skpt))
  
  # prior data conflict for enthusiastic prior
  box.enth <- sum(enth.post.nc[enth.post.nc <= enth.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(enth.post.nc)))
  print(paste0("Enthuastic prior compatibility: ", box.enth))
  
  # prior data conflict for non-informative prior
  box.ni <- sum(ni.post.nc[ni.post.nc <= ni.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with ni prior (should be 1): ", sum(ni.post.nc)))
  print(paste0("Non-informative prior compatibility: ", box.ni))
  
  # compute SKEPTICAL COMPONENT mixing weight
  eff.mix.prob <- 1 - max(box.enth - box.skpt, 0)
  return(cbind(eff.mix.prob, box.skpt, box.enth, box.ni))
}

start_time <- Sys.time()

x <- foreach (i = 1:nrow(Table0), .combine='c') %dopar% {
  prior_dat_conflict(Table0$y1.IP[i],
                     Table0$y0.IP[i],
                     Table0$y1.PC[i],
                     Table0$y0.PC[i])
}

x.t         <- data.frame(t(matrix(data = x, nrow = 4)))
names(x.t)  <- c("eff.mix.prob", "box.skpt", "box.enth", "box.ni")

end_time  <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
cat("Started  ", as.character(start_time), "\n",
    "Finished ", as.character(end_time), "\n",
    "Time difference of ", diff_time, " ", attr(diff_time, "units"), "\n",
    sep = "")

Table0 <- cbind(Table0, x.t)
Table0$risk.diff <- Table0$y1.IP/(Table0$y1.IP + Table0$y0.IP) - Table0$y1.PC/(Table0$y1.PC + Table0$y0.PC)
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################