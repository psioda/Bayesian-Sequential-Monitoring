ni_joint <- function(){
  
  q.outer    <- 1 - sig.eff
  scale      <- 1.5
  q.inner    <- (1 - pnorm(qnorm(1 - q.outer)/2) - q.outer)*scale
  
  source("priors/07_ni_marginal.R", local = TRUE)
  marginal()
  
  f1 <- function(a){   # q.outer
    ni.rd.prior <- function(x){
      dgnorm(x,           delta.ni.skpt, a[1], a[2])/
        (pgnorm(q = 1,    delta.ni.skpt, a[1], a[2]) -
           pgnorm(q = -1, delta.ni.skpt, a[1], a[2]))
    }
    
    # P(theta > theta_1) (called delta.ni.enth for historical reasons)
    integrate(ni.rd.prior, 
              lower = delta.ni.enth, 
              upper = 1)$value
  }
  
  f2 <- function(a){   # q.inner
    ni.rd.prior <- function(x){
      dgnorm(x,           delta.ni.skpt, a[1], a[2])/
        (pgnorm(q = 1,    delta.ni.skpt, a[1], a[2]) -
           pgnorm(q = -1, delta.ni.skpt, a[1], a[2]))
    }
    
    # P((theta_1+theta_0)/2 < theta < theta_1) (called delta.ni.intr for historical reasons)
    integrate(ni.rd.prior, 
              lower = delta.ni.intr, 
              upper = delta.ni.enth)$value
  }
  
  start <- c(ni.rd.alpha0,ni.rd.beta0)
  
  
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
  print(paste0("ni.rd.alpha0 ", a2[1]))
  print(paste0("ni.rd.beta0  ", a2[2]))
  
  ni.rd.alpha0 <- a2[1]
  ni.rd.beta0  <- a2[2]
  
  assign("ni.rd.alpha0",ni.rd.alpha0, envir = .GlobalEnv)
  assign("ni.rd.beta0", ni.rd.beta0,  envir = .GlobalEnv)
  
  source("priors/08_ni_conditional.R", local = TRUE)
  q.outer    <- 0.1
  scale      <- 1
  q.inner    <- (1 - pnorm(qnorm(1 - q.outer)/2) - q.outer)*scale
  mu1 <- mu + 0.2
  mu2 <- mu + 0.1
  ni_tail_area()
  
  f1 <- function(a){
    ni.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
      dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
        (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
           pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
        dgnorm(y,          mu, a[1], a[2])/
        (pgnorm(q = 1 - x, mu, a[1], a[2]) -
           pgnorm(q = 0,   mu, a[1], a[2]))
    }
    ni.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
      dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
        (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
           pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
        dgnorm(y,         mu,         a[1],    a[2])/
        (pgnorm(q = 1,    mu,         a[1],    a[2]) -
           pgnorm(q = -x, mu,         a[1],    a[2]))
    }
    # window around mu +/- delta.ni.enth
    c1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
    c2 <- integrate_debug(ni.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
    c3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
    c1 + c2 + c3
  }
  
  f2 <- function(a){
    ni.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
      dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
        (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
           pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
        dgnorm(y,          mu,        a[1],           a[2])/
        (pgnorm(q = 1 - x, mu,        a[1],           a[2]) -
           pgnorm(q = 0,   mu,        a[1],           a[2]))
    }
    ni.prior.2 <- function(x, y){ # for x < 0 (gamma < 0)
      dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
        (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
           pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
        dgnorm(y,         mu,         a[1],           a[2])/
        (pgnorm(q = 1,    mu,         a[1],           a[2]) -
           pgnorm(q = -x, mu,         a[1],           a[2]))
    }
    # window around mu +/- delta.ni.intr
    c1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
    c2 <- integrate_debug(ni.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
    c3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
    d1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
    d2 <- integrate_debug(ni.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
    d3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
    d1 + d2 + d3 - c1 - c2 - c3
  }
  
  start <- c(ni.alpha0, ni.beta0)
  
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
  print(paste0("ni.alpha0 ", a2[1]))
  print(paste0("ni.beta0  ", a2[2]))
  
  ni.alpha0 <- a2[1]
  ni.beta0  <- a2[2]
  
  assign("ni.alpha0", ni.alpha0, envir = .GlobalEnv)
  assign("ni.beta0", ni.beta0, envir = .GlobalEnv)
  
  # assemble final prior
  ni.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
    dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
      (pgnorm(q = 1,     delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
         pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
      dgnorm(y,          mu,         ni.alpha0,    ni.beta0)/
      (pgnorm(q = 1 - x, mu,         ni.alpha0,    ni.beta0) -
         pgnorm(q = 0,   mu,         ni.alpha0,    ni.beta0))
  }
  ni.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
    dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)/
      (pgnorm(q = 1,    delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
         pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0))*
      dgnorm(y,         mu,         ni.alpha0,    ni.beta0)/
      (pgnorm(q = 1,    mu,         ni.alpha0,    ni.beta0) -
         pgnorm(q = -x, mu,         ni.alpha0,    ni.beta0))
  }
  
  ## CHECK MARGINAL CONDITIONS
  print(paste0("ni marginal upper tail ", 
               integrate_debug(ni.prior.1, xmin = delta.ni.intr, xmax = delta.ni.enth, ymin = 0, ymax = function(x) 1 - x)))
  print(paste0("ni marginal upper half ",
               integrate_debug(ni.prior.1, xmin = delta.ni.enth, xmax = 1,          ymin = 0, ymax = function(x) 1 - x)))
  
  ## CHECK CONDITIONAL CONDITIONS
  c1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu1, ymin = mu1,            ymax = function(x) 1 - x)
  c2 <- integrate_debug(ni.prior.2, xmin = -mu1, xmax = 0,       ymin = mu1,            ymax = 1)
  c3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu1,    ymin = function(x) -x, ymax = 1)
  d1 <- integrate_debug(ni.prior.1, xmin = 0,    xmax = 1 - mu2, ymin = mu2,            ymax = function(x) 1 - x)
  d2 <- integrate_debug(ni.prior.2, xmin = -mu2, xmax = 0,       ymin = mu2,            ymax = 1)
  d3 <- integrate_debug(ni.prior.2, xmin = -1,   xmax = -mu2,    ymin = function(x) -x, ymax = 1)
  print(paste0("ni conditional upper tail ",
               c1 + c2 + c3))
  print(paste0("ni conditional upper half ",
               d1 + d2 + d3 - c1 - c2 - c3))
  print(paste0("ni joint integral ",
               integrate_debug(ni.prior.1,  xmin = 0,  xmax = 1 , ymin = 0,               ymax = function(x) 1 - x) +
                 integrate_debug(ni.prior.2,  xmin = -1, xmax = 0,  ymin = function(x) -x,  ymax = 1)))
}