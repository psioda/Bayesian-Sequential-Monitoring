rm(list=ls())

library(pracma) # numerical integration
library(gnorm)  # generalized normal distribution 

# trial data (IP response rate 60%, PC response rate 40%)
y1.IP <- 12
y0.IP <- 8
y1.PC <- 8
y0.PC <- 12

# risk difference prior parameters (derived in separate program)
delta.skpt     <- 0
skpt.rd.alpha0 <- 0.05776004
skpt.rd.beta0  <- 1.282048
delta.enth     <- 0.12
enth.rd.alpha0 <- 0.08658609
enth.rd.beta0  <- 2

# placebo prior parameters (derived in separate program)
mu          <- 0.39
skpt.alpha0 <- 0.2158015
skpt.beta0  <- 1.920246
enth.alpha0 <- 0.2190863
enth.beta0  <- 1.944127


# integrate function (progressively gets weaker if errors appear)
integrate_debug <- function(fun, xmin, xmax, ymin, ymax){
  tryCatch(
    tryCatch(
      tryCatch(
        integral2(fun, xmin, xmax, ymin, ymax)$Q,
        error = function(e) integral2(fun, xmin, xmax, ymin, ymax, singular = T)$Q),
    error = function(e) integral2(fun, xmin, xmax, ymin, ymax, abstol = 1E-6)$Q),
  error = function(e) integral2(fun, xmin, xmax, ymin, ymax, abstol = 1E-4)$Q)
}

# Function to compute prior predictive probability of the data under the skeptical
# and enthusiastic priors and compute prior data conflict quantity, which is then
# used to determine mixing weight in adaptive prior.
#
# Note: "x" corresponds to theta (risk difference) 
#       "y" corresponds to eta0 (placebo response rate)
#
#
prior_dat_conflict <- function(y1.IP, y0.IP, y1.PC, y0.PC){
  
  skpt.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  enth.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  
  for (i in 1:(y1.IP + y0.IP + 1)){
    print(paste0("Outer loop ", i, " of " , y1.IP + y0.IP + 1))
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
    }
  }
  
  # prior data conflict for skeptical prior
  skpt.psi <- sum(skpt.post.nc[skpt.post.nc <= skpt.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(skpt.post.nc)))
  print(paste0("Skeptical prior compatibility: ", skpt.psi))
  
  # prior data conflict for enthusiastic prior
  enth.psi <- sum(enth.post.nc[enth.post.nc <= enth.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(enth.post.nc)))
  print(paste0("Enthuastic prior compatibility: ", enth.psi))
  
  # compute SKEPTICAL COMPONENT mixing weight
  eff.mix.prob <- 1 - max(enth.psi - skpt.psi, 0)
  return(eff.mix.prob)
}

# test function
prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)