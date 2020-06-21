skpt_joint <- function(){

q.outer    <- 1 - sig.eff
scale      <- 0.75
q.inner    <- (1 - pnorm(qnorm(1 - q.outer)/2) - q.outer)*scale

source("priors/04_skpt_marginal.R", local = TRUE)
marginal()

f1 <- function(a){   # q.outer
  skpt.rd.prior <- function(x){
    dgnorm(x,           delta.skpt, a[1], a[2])/
      (pgnorm(q = 1,    delta.skpt, a[1], a[2]) -
         pgnorm(q = -1, delta.skpt, a[1], a[2]))
  }
  
  # P(theta > theta_1) (called delta.enth for historical reasons)
  integrate(skpt.rd.prior, 
            lower = delta.enth, 
            upper = 1)$value
}

f2 <- function(a){   # q.inner
  skpt.rd.prior <- function(x){
    dgnorm(x,           delta.skpt, a[1], a[2])/
      (pgnorm(q = 1,    delta.skpt, a[1], a[2]) -
         pgnorm(q = -1, delta.skpt, a[1], a[2]))
  }
  
  # P((theta_1+theta_0)/2 < theta < theta_1) (called delta.intr for historical reasons)
  integrate(skpt.rd.prior, 
            lower = delta.intr, 
            upper = delta.enth)$value
}

start <- c(skpt.rd.alpha0,skpt.rd.beta0)


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

print(paste0("Upper tail ", f1(a2)))
print(paste0("Upper half ", f2(a2)))
print(paste0("skpt.rd.alpha0 ", a2[1]))
print(paste0("skpt.rd.beta0  ", a2[2]))

skpt.rd.alpha0 <- a2[1]
skpt.rd.beta0  <- a2[2]

assign("skpt.rd.alpha0",skpt.rd.alpha0, envir = .GlobalEnv)
assign("skpt.rd.beta0", skpt.rd.beta0,  envir = .GlobalEnv)

source("priors/05_skpt_conditional.R", local = TRUE)
q.outer    <- 0.1
scale      <- 1
q.inner    <- (1 - pnorm(qnorm(1 - q.outer)/2) - q.outer)*scale
mu1 <- mu + 0.2
mu2 <- mu + 0.1
skpt_tail_area()

f1 <- function(a){
  skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
    dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
      (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
         pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
      dgnorm(y,          mu, a[1], a[2])/
      (pgnorm(q = 1 - x, mu, a[1], a[2]) -
         pgnorm(q = 0,   mu, a[1], a[2]))
  }
  skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
    dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
      (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
         pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
      dgnorm(y,         mu,         a[1],    a[2])/
      (pgnorm(q = 1,    mu,         a[1],    a[2]) -
         pgnorm(q = -x, mu,         a[1],    a[2]))
  }
  # window around mu +/- delta.enth
  c1 <- integrate_debug(skpt.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(skpt.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  c3 <- integrate_debug(skpt.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  c1 + c2 + c3
  }

f2 <- function(a){
  skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
    dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
      (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
         pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
      dgnorm(y,          mu,        a[1],           a[2])/
      (pgnorm(q = 1 - x, mu,        a[1],           a[2]) -
         pgnorm(q = 0,   mu,        a[1],           a[2]))
  }
  skpt.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
    dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
      (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
         pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
      dgnorm(y,         mu,         a[1],           a[2])/
      (pgnorm(q = 1,    mu,         a[1],           a[2]) -
         pgnorm(q = -x, mu,         a[1],           a[2]))
  }
  # window around mu +/- delta.intr
  c1 <- integrate_debug(skpt.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(skpt.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  c3 <- integrate_debug(skpt.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  d1 <- integrate_debug(skpt.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
  d2 <- integrate_debug(skpt.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
  d3 <- integrate_debug(skpt.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
  d1 + d2 + d3 - c1 - c2 - c3
  }

start <- c(skpt.alpha0, skpt.beta0)

fn <- function(a){
  (f1(a) - q.outer)^2 + (f2(a) - q.inner)^2
}

# nlm.fit <- nlminb(start     = start, 
#                   objective = fn, 
#                   lower     = c(0,0), 
#                   upper     = c(Inf, Inf))
optim.fit <- optim(par    = start,
             fn           = fn,
             lower        = c(0, 0),
             upper        = c(Inf, Inf),
             method       = "L-BFGS-B")
a2 <- optim.fit$par

print(paste0("Upper tail ", f1(a2)))
print(paste0("Upper half ", f2(a2)))
print(paste0("skpt.alpha0 ", a2[1]))
print(paste0("skpt.beta0  ", a2[2]))

skpt.alpha0 <- a2[1]
skpt.beta0  <- a2[2]

assign("skpt.alpha0", skpt.alpha0, envir = .GlobalEnv)
assign("skpt.beta0", skpt.beta0, envir = .GlobalEnv)

# assemble final prior
skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
    (pgnorm(q = 1,     delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
       pgnorm(q = -1,  delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
    dgnorm(y,          mu,         skpt.alpha0,    skpt.beta0)/
    (pgnorm(q = 1 - x, mu,         skpt.alpha0,    skpt.beta0) -
       pgnorm(q = 0,   mu,         skpt.alpha0,    skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
    (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
       pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
    dgnorm(y,         mu,         skpt.alpha0,    skpt.beta0)/
    (pgnorm(q = 1,    mu,         skpt.alpha0,    skpt.beta0) -
       pgnorm(q = -x, mu,         skpt.alpha0,    skpt.beta0))
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
}