##################################
### Evan Kwiatkowski, Feb 2020
##################################
# Objective: based on current data, create the following: 
# 1) scaled (un-normalized) posterior densities 
#    skpt.post.sc
#    enth.post.sc
# 2) scale factor
#    sc
# 3) scaled normalizing constant for posterior density
#    skpt.nc.sc
#    enth.nc.sc
# 4) scaled (un-normalized) marginal distributions
#    skpt.post.x.sc
#    skpt.post.y.sc
#    enth.post.x.sc
#    enth.post.y.sc
# 5) posterior mixing weights
#    fut.skpt.wt
#    fut.enth.wt
#    eff.skpt.wt
#    eff.enth.wt
#    inf.skpt.wt
#    inf.enth.wt
# Notes: Anything with "post" refers to a function

# SECTION 1: PRIORS (normalized)


skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  dgnorm(x,         delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
    (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
       pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
    dgnorm(y,          mu, skpt.alpha0, skpt.beta0)/
    (pgnorm(q = 1 - x, mu, skpt.alpha0, skpt.beta0) -
       pgnorm(q = 0,   mu, skpt.alpha0, skpt.beta0))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  dgnorm(x,           delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
    (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
       pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
    dgnorm(y,         mu,         skpt.alpha0,    skpt.beta0)/
    (pgnorm(q = 1,    mu,         skpt.alpha0,    skpt.beta0) -
       pgnorm(q = -x, mu,         skpt.alpha0,    skpt.beta0))
}

enth.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu,enth.alpha0, enth.beta0)/
    (pgnorm(q = 1 - x, mu, enth.alpha0, enth.beta0) -
       pgnorm(q = 0,   mu, enth.alpha0, enth.beta0))
}
enth.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  dgnorm(x,           delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu,enth.alpha0, enth.beta0)/
    (pgnorm(q = 1,     mu, enth.alpha0, enth.beta0) -
       pgnorm(q = -x,  mu, enth.alpha0, enth.beta0))
}

# SECTION 1.1 FIND MLEs (assuming nonzero cell counts in each)
PC.mle        <- y1.PC/sum(y0.PC, y1.PC)
IP.mle        <- y1.IP/sum(y0.IP, y1.IP)
risk.diff.mle <- IP.mle - PC.mle

# SECTION 2: PRIOR MIXING WEIGHTS (only if eff.mix.prob == NA)
box.skpt <- NA
box.enth <- NA
box.ni   <- NA

# if (is.na(eff.mix.prob)){ # commented out 06-17-21
#   if (IP.mle >= PC.mle){
#   skpt.lik <- skpt.prior.1(IP.mle - PC.mle, PC.mle)
#   enth.lik <- enth.prior.1(IP.mle - PC.mle, PC.mle)
#   } else {
#   skpt.lik <- skpt.prior.2(IP.mle - PC.mle, PC.mle)
#   enth.lik <- enth.prior.2(IP.mle - PC.mle, PC.mle)
#   }
#   eff.mix.prob <- 0.25 + 0.75*(skpt.lik/sum(skpt.lik, enth.lik)) # 3-19-20 update
# }

prior_dat_conflict <- function(y1.IP, y0.IP, y1.PC, y0.PC){
  
  # skpt.post.nc       <- matrix(NA, 
  #                              nrow = y1.IP + y0.IP + 1, 
  #                              ncol = y1.PC + y0.PC + 1)
  enth.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  # ni.post.nc         <- matrix(NA, 
  #                              nrow = y1.IP + y0.IP + 1, 
  #                              ncol = y1.PC + y0.PC + 1)
  
  print(paste0("Interim analysis ", sum(y1.IP,y0.IP,y1.PC,y0.PC)))
  
  for (i in 1:(y1.IP + y0.IP + 1)){
    for (j in 1:(y1.PC + y0.PC + 1)){
      # skpt.post.1   <- function(x, y){
      #   exp(
      #     dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
      #       dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) +
      #       dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) -
      #       log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
      #             pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
      #       dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) -
      #       log(pgnorm(q = 1 - x,  mu, skpt.alpha0, skpt.beta0) -
      #             pgnorm(q = 0,  mu, skpt.alpha0, skpt.beta0))
      #   )
      # }
      # skpt.post.2   <- function(x, y){
      #   exp(
      #     dbinom(i - 1,   y1.IP + y0.IP, x + y, log = TRUE) +
      #       dbinom(j - 1, y1.PC + y0.PC, y, log = TRUE) +
      #       dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) -
      #       log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
      #             pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
      #       dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) -
      #       log(pgnorm(q = 1,    mu, skpt.alpha0, skpt.beta0) -
      #             pgnorm(q = -x, mu, skpt.alpha0, skpt.beta0))
      #   )
      # }
      # skpt.post.nc[i, j] <- integrate_debug(skpt.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
      #   integrate_debug(skpt.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
      
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
      
      # ni.post.1 <- function(x, y){ # for x > 0 (theta > 0)
      #   exp(
      #     dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
      #       dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
      #   dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
      #     log(pgnorm(q = 1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
      #        pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
      #     dgnorm(y,          mu,            ni.alpha0,    ni.beta0, log = TRUE) -
      #     log(pgnorm(q = 1-x,mu,            ni.alpha0,    ni.beta0) -
      #        pgnorm(q = 0,   mu,            ni.alpha0,    ni.beta0))
      #   )
      # }
      # ni.post.2 <- function(x, y){ # for x < 0 (theta < 0)
      #   exp(
      #     dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
      #       dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
      #   dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
      #     log(pgnorm(q = 1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
      #        pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) + 
      #     dgnorm(y,         mu,            ni.alpha0,    ni.beta0, log = TRUE) -
      #     log(pgnorm(q = 1, mu,            ni.alpha0,    ni.beta0) -
      #        pgnorm(q = -x, mu,            ni.alpha0,    ni.beta0))
      #   )
      # }
      # ni.post.nc[i, j] <- integrate_debug(ni.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
      #                     integrate_debug(ni.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
      
    }
  }
  
  # # prior data conflict for skeptical prior
  # box.skpt <- sum(skpt.post.nc[skpt.post.nc <= skpt.post.nc[y1.IP + 1, y1.PC + 1]])
  # print(paste0("Sum marginal prob of data with skpt prior (should be 1): ", sum(skpt.post.nc)))
  # print(paste0("Skeptical prior compatibility: ", box.skpt))
  
  # prior data conflict for enthusiastic prior
  box.enth <- sum(enth.post.nc[enth.post.nc <= enth.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(enth.post.nc)))
  print(paste0("Enthuastic prior compatibility: ", box.enth))
  
  # # prior data conflict for non-informative prior
  # box.ni <- sum(ni.post.nc[ni.post.nc <= ni.post.nc[y1.IP + 1, y1.PC + 1]])
  # print(paste0("Sum marginal prob of data with ni prior (should be 1): ", sum(ni.post.nc)))
  # print(paste0("Non-informative prior compatibility: ", box.ni))
  
  # compute SKEPTICAL COMPONENT mixing weight
  # eff.mix.prob <- 1 - max(box.enth - box.skpt, 0)
  # eff.mix.prob20 <- max(box.skpt - box.enth, 0) # 6-10-21 edit
  
  # print(paste0("Efficacy mixing weight for skeptical component: ", eff.mix.prob))
  
  return(cbind(box.skpt, box.enth, box.ni))
  
  # return(cbind(eff.mix.prob, eff.mix.prob20, box.skpt, box.enth, box.ni))
}

delta.list <- seq(0, 0.25, by = 0.05)
beta.list <- c(1, 1.3220, 1.7370, 2.3220, 3.3220)
adapt.mat <- expand.grid(delta.list, beta.list)
adapt.mat$id <- 101:130

if (eff.mix.prob > 100){
  delta.a <- adapt.mat[adapt.mat$id == eff.mix.prob, "Var1"]
  beta.a <- adapt.mat[adapt.mat$id == eff.mix.prob, "Var2"]
  
  prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
  
  box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
  box.enth                   <- prior_data_conflict_result[,"box.enth"]
  box.ni                     <- prior_data_conflict_result[,"box.ni"]
  BoxPE                      <- box.enth
  eff.mix.prob               <- 1 - (1 - delta.a)*pbeta(BoxPE, 1, beta.a) # amount assigned to skeptical prior
}

# if (eff.mix.prob == 10){
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   eff.mix.prob               <- prior_data_conflict_result[,"eff.mix.prob"]
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
# }
# 
# if (eff.mix.prob == 20){ # 6-10-21 edit
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   eff.mix.prob               <- prior_data_conflict_result[,"eff.mix.prob20"] # 6-10-21 edit
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
# }
# 
# if (eff.mix.prob == 21){ # 6-17-21 edit
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   eff.mix.prob               <- max(prior_data_conflict_result[,"eff.mix.prob20"], 0.2) # 6-17-21 edit
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
# }
# 
# if (eff.mix.prob == 22){ # 6-17-21 edit
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   eff.mix.prob               <- max(prior_data_conflict_result[,"eff.mix.prob20"], 0.15) # 6-17-21 edit
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
# }
# 
# if (eff.mix.prob == 23){ # 6-17-21 edit
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   eff.mix.prob               <- max(prior_data_conflict_result[,"eff.mix.prob20"], 0.1) # 6-17-21 edit
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
# }
# 
# if (eff.mix.prob == 24){ # 6-17-21 edit
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   eff.mix.prob               <- max(prior_data_conflict_result[,"eff.mix.prob20"], 0.05) # 6-17-21 edit
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
# }
# 
# if (eff.mix.prob == 30){
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
#   eff.mix.prob               <- box.skpt / (box.skpt + box.enth)
# }


# SECTION 3: POSTERIOR DENSITIES
# log (un-normalized) posterior density
# Deemed unnecessary on 5/7/2020

# scale factor: average of estimated maximum of log (un-normalized) posterior densities
# Deemed unnecessary on 5/7/2020

# scaled (un-normalized) posterior density
# Edited on 5/7/2020
# Still should be on log scale

if (skpt.rd.alpha0 == 1E3 & skpt.alpha0 == 1E3){ 
  skpt.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE)
    )
  }
  skpt.post.sc.2 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE)
    )
  }
  enth.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE)
    )
  }
  enth.post.sc.2 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE)
    )
  }
} else if (skpt.alpha0 == 1E3){
  skpt.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) -
        log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
              pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))
    )
  }
  skpt.post.sc.2 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) +
        dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) -
        log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
              pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))
    )
  }
  enth.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
              pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))
    )
  }
  enth.post.sc.2 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
              pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))
    )
  }
} else {
  skpt.post.sc.1 <- function(x, y){
    exp(
      dbinom(y1.IP,   y1.IP + y0.IP, x + y, log = TRUE) +
        dbinom(y1.PC, y1.PC + y0.PC, y, log = TRUE) + 
        dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
        log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
              pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
        dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
        log(pgnorm(q = 1-x,  mu, skpt.alpha0, skpt.beta0) - 
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
        log(pgnorm(q = 1-x,  mu, enth.alpha0, enth.beta0) - 
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
}
# scaled normalizing constant for posterior density
# No longer scaled on 5/7/2020, just regular probability of data
skpt.nc.sc <- integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
  integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
enth.nc.sc <- integrate_debug(enth.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
  integrate_debug(enth.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

# check that posterior is normalized
# print(paste0("skpt posterior integral: ",
#      (integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#       integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
#       skpt.nc.sc))
# print(paste0("enth posterior integral: ",
#      (integrate_debug(enth.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#       integrate_debug(enth.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1))/
#       enth.nc.sc))