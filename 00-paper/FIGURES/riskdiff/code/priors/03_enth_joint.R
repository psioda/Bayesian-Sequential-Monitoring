enth_joint <- function(){
  
q.outer    <- 1 - sig.fut
scale      <- 1
q.inner    <- (1 - pnorm(qnorm(1 - q.outer)/2) - q.outer)*scale 

source("priors/01_enth_marginal.R", local = TRUE)
marginal()

f1 <- function(a){   # q.outer
  enth.rd.prior <- function(x){
    dgnorm(x,           delta.enth, a[1], a[2])/
      (pgnorm(q = 1,    delta.enth, a[1], a[2]) -
         pgnorm(q = -1, delta.enth, a[1], a[2]))
  }
  
  # P(theta < theta_0) (called delta.skpt for historical reasons)
  integrate(enth.rd.prior, 
            lower = -1,         
            upper = delta.skpt)$value
}
f2 <- function(a){   # q.inner
  enth.rd.prior <- function(x){
    dgnorm(x,           delta.enth, a[1], a[2])/
      (pgnorm(q = 1,    delta.enth, a[1], a[2]) -
         pgnorm(q = -1, delta.enth, a[1], a[2]))
  }
  
  # P(theta_0 < theta < (theta_1+theta_0)/2) (called delta.intr for historical reasons)
  integrate(enth.rd.prior, 
            lower = delta.skpt, 
            upper = delta.intr)$value
}

start <- c(enth.rd.alpha0,enth.rd.beta0)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

optim.fit <- optim(par          = start,
                   fn           = fn,
                   lower        = c(0, 0),
                   upper        = c(Inf, Inf),
                   method       = "L-BFGS-B")

a2 <- optim.fit$par
nlm.fit <- nlminb(start     = start,
                  objective = fn,
                  lower     = c(0, 0),
                  upper     = c(Inf, Inf))
a2 <- nlm.fit$par

print(paste0("Lower tail ", f1(a2)))
print(paste0("Lower half ", f2(a2)))
print(paste0("enth.rd.alpha0 ", a2[1]))
print(paste0("enth.rd.beta0  ", a2[2]))

enth.rd.alpha0 <- a2[1]
enth.rd.beta0  <- a2[2]

assign("enth.rd.alpha0", enth.rd.alpha0, envir = .GlobalEnv)
assign("enth.rd.beta0",  enth.rd.beta0,  envir = .GlobalEnv)

source("priors/02_enth_conditional.R", local = TRUE)
q.outer    <- 0.1
scale      <- 1
q.inner    <- (1 - pnorm(qnorm(1 - q.outer)/2) - q.outer)*scale
mu1 <- mu + 0.2
mu2 <- mu + 0.1
enth_tail_area()

f1 <- function(a){
  enth.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
    dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
      (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
         pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
      dgnorm(y,          mu, a[1], a[2])/
      (pgnorm(q = 1 - x, mu, a[1], a[2]) -
         pgnorm(q = 0,   mu, a[1], a[2]))
  }
  enth.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
    dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
      (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
         pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
      dgnorm(y,         mu, a[1], a[2])/
      (pgnorm(q = 1,    mu, a[1], a[2]) -
         pgnorm(q = -x, mu, a[1], a[2]))
  }
  # window around mu +/- delta.enth
  c1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(enth.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  
  if (enth.prior.2(mu1,-mu1) > 1E-10){
  c3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  } else {
  c3 <- 0
  }
  c1 + c2 + c3
  }

f2 <- function(a){
  enth.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
    dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
      (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
         pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
      dgnorm(y,          mu, a[1], a[2])/
      (pgnorm(q = 1 - x, mu, a[1], a[2]) -
         pgnorm(q = 0,   mu, a[1], a[2]))
  }
  enth.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
    dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
      (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
         pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
      dgnorm(y,         mu, a[1], a[2])/
      (pgnorm(q = 1,    mu, a[1], a[2]) -
         pgnorm(q = -x, mu, a[1], a[2]))
  }
  # window around mu +/- delta.intr
  c1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(enth.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  if (enth.prior.2(-mu1, mu1) > 1E-10){
    c3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  } else {
    c3 <- 0
  }
  d1 <- integrate_debug(enth.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
  d2 <- integrate_debug(enth.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
  if (enth.prior.2(-mu2, mu2) > 1E-10){
    d3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
  } else {
    d3 <- 0
  }
  d1 + d2 + d3 - c1 - c2 - c3
  }

start <- c(enth.alpha0, enth.beta0)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

nlm.fit <- nlminb(start     = start, 
                   objective = fn, 
                   lower     = c(0,0), 
                   upper     = c(Inf, Inf))

optim.fit <- optim(par    = start,
             fn           = fn,
             lower        = c(0, 0),
             upper        = c(Inf, Inf),
             method       = "L-BFGS-B")
a2 <- optim.fit$par

print(paste0("Upper tail ", f1(a2)))
print(paste0("Upper half ", f2(a2)))
print(paste0("enth.alpha0 ", a2[1]))
print(paste0("enth.beta0  ", a2[2]))

enth.alpha0 <- a2[1]
enth.beta0  <- a2[2]

assign("enth.alpha0", enth.alpha0, envir = .GlobalEnv)
assign("enth.beta0",  enth.beta0,  envir = .GlobalEnv)

# assemble final prior
enth.prior.1 <- function(x, y){ # for x > 0 (gamma > 0)
  dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu, enth.alpha0, enth.beta0)/
    (pgnorm(q = 1 - x, mu, enth.alpha0, enth.beta0) -
       pgnorm(q = 0,   mu, enth.alpha0, enth.beta0))
}
enth.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
  dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu, enth.alpha0, enth.beta0)/
    (pgnorm(q = 1,     mu, enth.alpha0, enth.beta0) -
       pgnorm(q = -x,  mu, enth.alpha0, enth.beta0))
}

## CHECK MARGINAL CONDITIONS
print(paste0("Enth marginal lower tail ", 
             integrate_debug(enth.prior.1, 
                             xmin = -1,         
                             xmax = delta.skpt, # careful, only one piece since delta.skpt = 0
                             ymin = 0, 
                             ymax = function(x) 1 - x)))
print(paste0("Enth marginal lower half ",
             integrate_debug(enth.prior.1, 
                             xmin = delta.skpt, 
                             xmax = delta.intr, 
                             ymin = 0, 
                             ymax = function(x) 1 - x)))

## CHECK CONDITIONAL CONDITIONS
c1 <- integrate_debug(enth.prior.1,   xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
c2 <- integrate_debug(enth.prior.2,   xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
if (enth.prior.2(-mu1, mu1) > 1E-10){
  c3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
} else {
  c3 <- 0
}
d1 <- integrate_debug(enth.prior.1,   xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
d2 <- integrate_debug(enth.prior.2,   xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
if (enth.prior.2(-mu2, mu2) > 1E-10){
  d3 <- integrate_debug(enth.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
} else {
  d3 <- 0
}
print(paste0("Enth conditional upper tail ",
             c1 + c2 + c3))
print(paste0("Enth conditional upper half ",
             d1 + d2 + d3 - c1 - c2 - c3))
print(paste0("Enth joint integral ",
             integrate_debug(enth.prior.1,  xmin = 0,  xmax = 1 , ymin = 0,               ymax = function(x) 1 - x) +
               integrate_debug(enth.prior.2,  xmin = -1, xmax = 0,  ymin = function(x) -x,  ymax = 1)))
}

