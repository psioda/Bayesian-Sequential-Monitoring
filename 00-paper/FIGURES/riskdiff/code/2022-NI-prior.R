rm(list = ls())
load(file = 'args_model.RData') # loads all model information include prior parameters AND SETS SEED
# skpt.alpha0    <- 1E3     # eta
skpt.rd.alpha0 <- 100  # theta
mu1 <- mu + 0.2
mu2 <- mu + 0.1

# mu1 <- 0.8
# mu2 <- 0
# mu <- 0.8

# # # December 2021

# skpt.rd.alpha0 <- 100  # theta
# skpt.alpha0 <- 100
# skpt.rd.alpha0 <- 1000

# # assemble final prior
# skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
#   dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
#     (pgnorm(q = 1,     delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
#        pgnorm(q = -1,  delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
#     dgnorm(y,          mu,         skpt.alpha0,    skpt.beta0)/
#     (pgnorm(q = 1 - x, mu,         skpt.alpha0,    skpt.beta0) -
#        pgnorm(q = 0,   mu,         skpt.alpha0,    skpt.beta0))
# }
# skpt.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
#   dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
#     (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
#        pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
#     dgnorm(y,         mu,         skpt.alpha0,    skpt.beta0)/
#     (pgnorm(q = 1,    mu,         skpt.alpha0,    skpt.beta0) -
#        pgnorm(q = -x, mu,         skpt.alpha0,    skpt.beta0))
# }

skpt.prior.1 <- function(x, y){
  x*0 + y*0 + 1
}

skpt.prior.2 <- function(x, y){
  x*0 + y*0 + 1
}

## CHECK MARGINAL CONDITIONS
print(paste0("Skpt marginal upper tail ", 
             integrate_debug(skpt.prior.1, xmin = delta.intr, xmax = delta.enth, ymin = 0, ymax = function(x) 1 - x)))
print(paste0("Skpt marginal upper half ",
             integrate_debug(skpt.prior.1, xmin = delta.enth, xmax = 1,          ymin = 0, ymax = function(x) 1 - x)))

## CHECK CONDITIONAL CONDITIONS
c1 <- integrate_debug(skpt.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
c2 <- integrate_debug(skpt.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
c3 <- integrate_debug(skpt.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
d1 <- integrate_debug(skpt.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
d2 <- integrate_debug(skpt.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
d3 <- integrate_debug(skpt.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
print(paste0("Skpt conditional upper tail ",
             c1 + c2 + c3))
print(paste0("Skpt conditional upper half ",
             d1 + d2 + d3 - c1 - c2 - c3))
print(paste0("Skpt joint integral ",
             integrate_debug(skpt.prior.1,  xmin = 0,  xmax = 1 , ymin = 0,               ymax = function(x) 1 - x) +
               integrate_debug(skpt.prior.2,  xmin = -1, xmax = 0,  ymin = function(x) -x,  ymax = 1)))

