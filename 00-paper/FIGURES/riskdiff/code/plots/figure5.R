#########################################
#### Figure 5, Risk Diff Prior Plots ####
#########################################

## SKEPTICAL PRIORS ##
# (un-normalized) prior density
skpt.prior <- function(x, y){
  exp(
    - (abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
      (abs((y - x) - delta.skpt)/skpt.sigma0)^skpt.lambda0
  )
}
# normalizing constant for prior density
skpt.prior.nc <- integrate_debug(fun = skpt.prior,
                                 xmin = 0,
                                 xmax = 1,
                                 ymin = 0,
                                 ymax = 1)

## ENTHUASTIC PRIORS ##
# (un-normalized) prior density
enth.prior <- function(x, y){
  exp(
    - (abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
      (abs((y - x) - delta.enth)/enth.sigma0)^enth.lambda0
  )
}

# normalizing constant for prior density
enth.prior.nc <- integrate_debug(fun = enth.prior,
                                 xmin = 0,
                                 xmax = 1,
                                 ymin = 0,
                                 ymax = 1)

x.len <- 101
grid1d <- seq(0, 1, length = x.len)
x <- grid1d
y <- grid1d
grid <- expand.grid(x = x, y = y)

grid.skpt <- grid
grid.skpt$z <- skpt.prior(grid.skpt$x, grid.skpt$y)/skpt.prior.nc

grid.enth <- grid
grid.enth$z <- enth.prior(grid.enth$x, grid.enth$y)/enth.prior.nc

# need to get distribution of theta1 - theta0 (y - x)
# uses convolution formula https://www.math.arizona.edu/~jwatkins/n-bivariate.pdf (typos in example 3)
# u = x - y, v = x

# (un-normalized) marginal density
marginal.u.skpt <- function(v){
  exp(
    -(abs(v - mu)/placebo.sigma0)^placebo.lambda0 -
      (abs(u - delta.skpt)/skpt.sigma0)^skpt.lambda0
  )
}

# (un-normalized) marginal density
marginal.u.enth <- function(v){
  exp(
    -(abs(v - mu)/placebo.sigma0)^placebo.lambda0 -
      (abs(u - delta.enth)/enth.sigma0)^enth.lambda0
  )
}

# rd = risk difference

rd.length <- 101
rd.grid.lower <- seq(-1, 0, length = rd.length)
rd.grid.upper <- seq(0,  1, length = rd.length)

skpt.lower <- NA
enth.lower <- NA
skpt.upper <- NA
enth.upper <- NA

for (i in 1:rd.length){
  u <- rd.grid.lower[i]
  skpt.lower[i] <- (integrate(marginal.u.skpt, lower = -u, upper = 1)$value)/skpt.prior.nc # fixed typo 1/30/20
  enth.lower[i] <- (integrate(marginal.u.enth, lower = -u, upper = 1)$value)/enth.prior.nc # fixed typo 1/30/20
  u <- rd.grid.upper[i]
  skpt.upper[i] <- (integrate(marginal.u.skpt, lower = 0, upper = 1 - u)$value)/skpt.prior.nc # fixed typo 1/30/20
  enth.upper[i] <- (integrate(marginal.u.enth, lower = 0, upper = 1 - u)$value)/enth.prior.nc # fixed typo 1/30/20
}

# check integrals of desired regions
trapz(c(rd.grid.lower, rd.grid.upper),
      c(skpt.lower, skpt.upper))
trapz(c(rd.grid.upper[rd.grid.upper >= delta.enth]),
      c(skpt.upper[rd.grid.upper >= delta.enth]))

trapz(c(rd.grid.lower,rd.grid.upper), 
      c(enth.lower,enth.upper))
trapz(c(rd.grid.lower),
      c(enth.lower))

skpt.x<-0
skpt.y<-0
enth.x<-0
enth.y<-0

for (i in 1:length(grid1d)){
  ## SKEPTICAL PRIORS ##
  x<-grid1d[i]
  marginal.x <- function(y){
    exp(
      -(abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
      (abs((y-x) - delta.skpt)/skpt.sigma0)^skpt.lambda0
    )
    }
  skpt.x[i]<-(integrate(marginal.x,lower=0,upper=1)$value)/skpt.prior.nc
  
  y<-grid1d[i]
  marginal.y <- function(x){
    exp(-(abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
          (abs((y - x)-delta.skpt)/skpt.sigma0)^skpt.lambda0
    )
    }
  skpt.y[i] <- (integrate(marginal.y,lower=0,upper=1)$value)/skpt.prior.nc
  
  ## ENTHUASTIC PRIORS ##
  x <- grid1d[i]
  marginal.x<-function(y){
    exp(
      -(abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
        (abs((y - x)-delta.enth)/enth.sigma0)^enth.lambda0
    )
    }
  enth.x[i] <- (integrate(marginal.x,lower=0,upper=1)$value)/enth.prior.nc
  
  y <- grid1d[i]
  marginal.y <- function(x){
    exp(-(abs(x - mu)/placebo.sigma0)^placebo.lambda0 - 
          (abs((y - x)-delta.enth)/enth.sigma0)^enth.lambda0
    )
    }
  enth.y[i] <- (integrate(marginal.y,lower=0,upper=1)$value)/enth.prior.nc
} 