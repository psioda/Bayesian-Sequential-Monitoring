library(foreach)
library(doParallel)
# register cores from doParallel package
registerDoParallel(detectCores())
getDoParWorkers()

setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/Real FDA Data Example/MP_FDA_Check/output-inference/")

# foreach(y1.IP = 29:59) %dopar% {

for(y1.IP in 15:30){
  
  # y1.IP <- 45
  y0.IP <- (50 - y1.IP)
  y1.PC <- 20
  y0.PC <- 30
  
  # ss <- 1E2
  # PC.p <- 0.39
  # IP.p <- 0.45
  # y1.IP <- IP.p*ss
  # y0.IP <- ss - y1.IP
  # y1.PC <- PC.p*ss
  # y0.PC <- ss - y1.PC
  
  # y1.IP <- dat[dat$targOutNum == index, "yObs1"]
  # y0.IP <- dat[dat$targOutNum == index, "nObs1"] - dat[dat$targOutNum == index, "yObs1"]
  # y1.PC <- dat[dat$targOutNum == index, "yObs0"]
  # y0.PC <- dat[dat$targOutNum == index, "nObs0"] - dat[dat$targOutNum == index, "yObs0"]
  
  # y1.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 1)
  # y0.IP <- sum(responses.IP[outcome.times.IP <= outcome.times.all[index]] == 0)
  # y1.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 1)
  # y0.PC <- sum(responses.PC[outcome.times.PC <= outcome.times.all[index]] == 0)
  
  # SECTION 1.1 FIND MLEs (assuming nonzero cell counts in each)
  PC.mle        <- y1.PC/sum(y0.PC, y1.PC)
  IP.mle        <- y1.IP/sum(y0.IP, y1.IP)
  risk.diff.mle <- IP.mle - PC.mle
  
  # SECTION 2: PRIOR MIXING WEIGHTS (only if eff.mix.prob == NA)
  box.skpt <- NA
  box.enth <- NA
  box.ni   <- NA
  
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
    
    print(paste0("Interim analysis ", sum(y1.IP,y0.IP,y1.PC,y0.PC)))
    
    for (i in 1:(y1.IP + y0.IP + 1)){
      for (j in 1:(y1.PC + y0.PC + 1)){
        skpt.post.1   <- function(x, y){
          exp(
            dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) +
              dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) -
              log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
                    pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
              dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) -
              log(pgnorm(q = 1 - x,  mu, skpt.alpha0, skpt.beta0) -
                    pgnorm(q = 0,  mu, skpt.alpha0, skpt.beta0))
          )
        }
        skpt.post.2   <- function(x, y){
          exp(
            dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) +
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
            dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) + 
              dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
              log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
                    pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
              dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
              log(pgnorm(q = 1 - x,  mu, enth.alpha0, enth.beta0) - 
                    pgnorm(q = 0,  mu, enth.alpha0, enth.beta0))
          )
        }
        enth.post.2   <- function(x, y){
          exp(
            dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) + 
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
            dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) +
              dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
              log(pgnorm(q = 1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
                    pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
              dgnorm(y,          mu,            ni.alpha0,    ni.beta0, log = TRUE) -
              log(pgnorm(q = 1 - x,mu,            ni.alpha0,    ni.beta0) -
                    pgnorm(q = 0,   mu,            ni.alpha0,    ni.beta0))
          )
        }
        ni.post.2 <- function(x, y){ # for x < 0 (theta < 0)
          exp(
            dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) +
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
    
    # # prior data conflict for skeptical prior
    box.skpt <- sum(skpt.post.nc[skpt.post.nc <= skpt.post.nc[y1.IP + 1, y1.PC + 1]])
    print(paste0("Sum marginal prob of data with skpt prior (should be 1): ", sum(skpt.post.nc)))
    print(paste0("Skeptical prior compatibility: ", box.skpt))
    
    # prior data conflict for enthusiastic prior
    box.enth <- sum(enth.post.nc[enth.post.nc <= enth.post.nc[y1.IP + 1, y1.PC + 1]])
    print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(enth.post.nc)))
    print(paste0("Enthuastic prior compatibility: ", box.enth))
    
    # # prior data conflict for non-informative prior
    box.ni <- sum(ni.post.nc[ni.post.nc <= ni.post.nc[y1.IP + 1, y1.PC + 1]])
    print(paste0("Sum marginal prob of data with ni prior (should be 1): ", sum(ni.post.nc)))
    print(paste0("Non-informative prior compatibility: ", box.ni))
    
    return(cbind(box.skpt, box.enth, box.ni))
  }
  
  prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
  box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
  box.enth                   <- prior_data_conflict_result[,"box.enth"]
  box.ni                     <- prior_data_conflict_result[,"box.ni"]
  
  # SECTION 3: POSTERIOR DENSITIES
  skpt.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
              pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
        dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
        log(pgnorm(q = 1 - x,  mu, skpt.alpha0, skpt.beta0) - 
              pgnorm(q = 0,  mu, skpt.alpha0, skpt.beta0))
    )
  }
  skpt.post.sc.2 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
              pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
        dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    mu, skpt.alpha0, skpt.beta0) - 
              pgnorm(q = -x, mu, skpt.alpha0, skpt.beta0))
    )
  }
  enth.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
              pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
        dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
        log(pgnorm(q = 1 - x,  mu, enth.alpha0, enth.beta0) - 
              pgnorm(q = 0,  mu, enth.alpha0, enth.beta0))
    )
  }
  enth.post.sc.2 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
              pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
        dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    mu, enth.alpha0, enth.beta0) - 
              pgnorm(q = -x, mu, enth.alpha0, enth.beta0))
    )
  }
  ni.post.sc.1 <- function(x, y){ # for x > 0 (theta > 0)
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
        log(pgnorm(q = 1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
              pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
        dgnorm(y,          mu,            ni.alpha0,    ni.beta0, log = TRUE) -
        log(pgnorm(q = 1 - x,mu,            ni.alpha0,    ni.beta0) -
              pgnorm(q = 0,   mu,            ni.alpha0,    ni.beta0))
    )
  }
  ni.post.sc.2 <- function(x, y){ # for x < 0 (theta < 0)
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
        log(pgnorm(q = 1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
              pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
        dgnorm(y,         mu,            ni.alpha0,    ni.beta0, log = TRUE) -
        log(pgnorm(q = 1, mu,            ni.alpha0,    ni.beta0) -
              pgnorm(q = -x, mu,            ni.alpha0,    ni.beta0))
    )
  }
  
  # scaled normalizing constant for posterior density
  # No longer scaled on 5/7/2020, just regular probability of data
  skpt.nc.sc <- integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
    integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
  enth.nc.sc <- integrate_debug(enth.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
    integrate_debug(enth.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
  ni.nc.sc <- integrate_debug(ni.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
    integrate_debug(ni.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
  
  # POSTERIOR DENSITIES FOR EXPECTED VALUE
  skpt.post.sc.1.x <- function(x, y){
    x * exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
              pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
        dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
        log(pgnorm(q = 1 - x,  mu, skpt.alpha0, skpt.beta0) - 
              pgnorm(q = 0,  mu, skpt.alpha0, skpt.beta0))
    )
  }
  skpt.post.sc.2.x <- function(x, y){
    x * exp(dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
              dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
              dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
              log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
                    pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
              dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
              log(pgnorm(q = 1,    mu, skpt.alpha0, skpt.beta0) - 
                    pgnorm(q = -x, mu, skpt.alpha0, skpt.beta0))
    )
  }
  enth.post.sc.1.x <- function(x, y){
    x * exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
              pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
        dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
        log(pgnorm(q = 1 - x,  mu, enth.alpha0, enth.beta0) - 
              pgnorm(q = 0,  mu, enth.alpha0, enth.beta0))
    )
  }
  enth.post.sc.2.x <- function(x, y){
    x * exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
              pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
        dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    mu, enth.alpha0, enth.beta0) - 
              pgnorm(q = -x, mu, enth.alpha0, enth.beta0))
    )
  }
  ni.post.sc.1.x <- function(x, y){ # for x > 0 (theta > 0)
    x * exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
        log(pgnorm(q = 1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
              pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
        dgnorm(y,          mu,            ni.alpha0,    ni.beta0, log = TRUE) -
        log(pgnorm(q = 1 - x,mu,            ni.alpha0,    ni.beta0) -
              pgnorm(q = 0,   mu,            ni.alpha0,    ni.beta0))
    )
  }
  ni.post.sc.2.x <- function(x, y){ # for x < 0 (theta < 0)
    x * exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
        log(pgnorm(q = 1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
              pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
        dgnorm(y,         mu,            ni.alpha0,    ni.beta0, log = TRUE) -
        log(pgnorm(q = 1, mu,            ni.alpha0,    ni.beta0) -
              pgnorm(q = -x, mu,            ni.alpha0,    ni.beta0))
    )
  }
  
  # Compute posterior means
  skpt.pm.x <- (integrate_debug(skpt.post.sc.1.x, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                  integrate_debug(skpt.post.sc.2.x, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)) / skpt.nc.sc
  enth.pm.x <- (integrate_debug(enth.post.sc.1.x, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                  integrate_debug(enth.post.sc.2.x, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)) / enth.nc.sc
  ni.pm.x <- (integrate_debug(ni.post.sc.1.x, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
                integrate_debug(ni.post.sc.2.x, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)) / ni.nc.sc
  
  # Single inference priors
  pm.skpt <- skpt.pm.x
  pm.enth <- enth.pm.x
  pm.ni <- ni.pm.x
  pm.rd <- risk.diff.mle
  
  # Fixed weight inference prior
  omega <- 0.5*skpt.nc.sc/(0.5*skpt.nc.sc + 0.5*enth.nc.sc)
  pm.50.50 <- omega*skpt.pm.x + (1 - omega)*enth.pm.x
  
  # 2-part adaptive weight inference prior
  box.weighted.skpt <- box.skpt / (box.skpt + box.enth)
  omega <- box.weighted.skpt*skpt.nc.sc / (box.weighted.skpt*skpt.nc.sc + (1 - box.weighted.skpt)*enth.nc.sc)
  pm.skpt.enth <- omega*skpt.pm.x + (1 - omega)*enth.pm.x
  
  # 3-part adaptive weight inference prior
  box.weighted.skpt <- box.skpt / (box.skpt + box.enth + box.ni)
  box.weighted.enth <- box.enth / (box.skpt + box.enth + box.ni)
  box.weighted.ni   <- box.ni / (box.skpt + box.enth + box.ni)
  omega.skpt <- box.weighted.skpt*skpt.nc.sc / (box.weighted.skpt*skpt.nc.sc + box.weighted.enth*enth.nc.sc + box.weighted.ni*ni.nc.sc)
  omega.enth <- box.weighted.enth*enth.nc.sc / (box.weighted.skpt*skpt.nc.sc + box.weighted.enth*enth.nc.sc + box.weighted.ni*ni.nc.sc)
  omega.ni   <- box.weighted.ni*ni.nc.sc     / (box.weighted.skpt*skpt.nc.sc + box.weighted.enth*enth.nc.sc + box.weighted.ni*ni.nc.sc)
  pm.skpt.enth.ni   <- omega.skpt*skpt.pm.x + omega.enth*enth.pm.x + omega.ni*ni.pm.x
  
  # 3-part adaptive weight inference prior w/ excess probability
  box.ni.tilde <- max(0, box.ni - max(box.skpt,box.enth))
  box.weighted.skpt <- box.skpt / (box.skpt + box.enth + box.ni.tilde)
  box.weighted.enth <- box.enth / (box.skpt + box.enth + box.ni.tilde)
  box.weighted.ni   <- box.ni.tilde / (box.skpt + box.enth + box.ni.tilde)
  omega.skpt <- box.weighted.skpt*skpt.nc.sc / (box.weighted.skpt*skpt.nc.sc + box.weighted.enth*enth.nc.sc + box.weighted.ni*ni.nc.sc)
  omega.enth <- box.weighted.enth*enth.nc.sc / (box.weighted.skpt*skpt.nc.sc + box.weighted.enth*enth.nc.sc + box.weighted.ni*ni.nc.sc)
  omega.ni   <- box.weighted.ni*ni.nc.sc     / (box.weighted.skpt*skpt.nc.sc + box.weighted.enth*enth.nc.sc + box.weighted.ni*ni.nc.sc)
  pm.skpt.enth.ni.tilde <- omega.skpt*skpt.pm.x + omega.enth*enth.pm.x + omega.ni*ni.pm.x
  
  write.csv(data.frame(cbind(pm.skpt,
                             pm.enth,
                             pm.ni,
                             pm.rd,
                             pm.50.50,
                             pm.skpt.enth,
                             pm.skpt.enth.ni,
                             pm.skpt.enth.ni.tilde,
                             box.skpt, box.enth, box.ni, 
                             y1.IP, y1.PC, y0.IP, y0.PC)), file = paste0(y1.IP, ".csv"))
  
  # return(data.frame(cbind(pm.skpt,
  #                         pm.enth,
  #                         pm.ni,
  #                         pm.rd,
  #                         pm.50.50,
  #                         pm.skpt.enth,
  #                         pm.skpt.enth.ni,
  #                         pm.skpt.enth.ni.tilde,
  #                         box.skpt, box.enth, box.ni, 
  #                         y1.IP, y1.PC, y0.IP, y0.PC)))
  
  
}