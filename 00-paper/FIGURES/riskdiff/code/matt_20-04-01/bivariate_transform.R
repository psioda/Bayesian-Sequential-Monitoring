rm(list = ls())
load(file = 'P:/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/args_model.RData')
library(gnorm)
library(pracma)

y0.PC <- 10
y0.IP <- 10
y1.PC <- 10
y1.IP <- 10

PC.mle <- y1.PC/sum(y0.PC, y1.PC)
IP.mle <- y1.IP/sum(y0.IP, y1.IP)

# SECTION 3: POSTERIOR DENSITIES
# log (un-normalized) posterior density
skpt.post.log.1 <- function(x, y){
  y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
  y1.PC*log(y)     + y0.PC*log(1 - y) -
    (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
    (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
    log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}
skpt.post.log.2 <- function(x, y){
  y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
  y1.PC*log(y)     + y0.PC*log(1 - y) -
    (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
    (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
    log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
}

# scale factor: average of estimated maximum of log (un-normalized) posterior densities
if (IP.mle >= PC.mle){
  sc <- (skpt.post.log.1(IP.mle - PC.mle, PC.mle))
} else {
  sc <- (skpt.post.log.2(IP.mle - PC.mle, PC.mle))
}

# scaled (un-normalized) posterior density
skpt.post.sc.1 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}
skpt.post.sc.2 <- function(x, y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}

skpt.post.sc.12 <- function(x,y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
    y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1*(1-x > 1) + (1-x)*(1-x <= 1), mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -x*(x<=0),                      mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}

integrate_debug(skpt.post.sc.12, xmin = 0, xmax = 1, ymin = 0, ymax = function(x) 1 - x)
integrate_debug(skpt.post.sc.1,  xmin = 0, xmax = 1, ymin = 0, ymax = function(x) 1 - x)

integrate_debug(skpt.post.sc.12, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
integrate_debug(skpt.post.sc.2,  xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1) 

skpt.post.sc.12a <- function(x,y){
  exp(
    y1.IP*log(x) + y0.IP*log(1 - (x)) + 
      y1.PC*log(y)     + y0.PC*log(1 - y) -
      (abs((x - y) - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1*(1-(x - y) > 1) + (1-(x - y))*(1-(x - y) <= 1), mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -(x - y)*((x - y)<=0),                            mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}

integrate_debug(skpt.post.sc.12a, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
skpt.post.sc.12a(0.6,0.6)




bivariate.test <- function(x,y){
  exp(
    y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
    y1.PC*log(x)     + y0.PC*log(1 - x) -
      (abs(y - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(x - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1*(1-y > 1) + (1 - y)*(1 - y <= 1), mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -y*(y<=0),                     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}

integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = function(x) 1 - x)

integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = 0)
integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = 0,              ymax = function(x) 1 - x)
integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x)

integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = function(x) 1 - x)
integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = 0) + 
integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = 0,              ymax = function(x) 1 - x)

bivariate.test2 <- function(x,y){
  exp(
    y1.IP*log(y) + y0.IP*log(1 - (y)) + 
      y1.PC*log(x)     + y0.PC*log(1 - x) -
      (abs(y - x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
      (abs(x - mu)/skpt.alpha0)^skpt.beta0 -
      log(pgnorm(q = 1*(1-(y - x) > 1) + (1 - (y - x))*(1 - (y - x) <= 1), mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
            pgnorm(q = -(y - x)*((y - x)<=0),                             mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
      sc
  )
}

integrate_debug(bivariate.test2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
integrate_debug(bivariate.test2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = 0) + 
  integrate_debug(bivariate.test, xmin = 0, xmax = 1, ymin = 0,              ymax = function(x) 1 - x)


integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x)
integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)


# scaled normalizing constant for posterior density
skpt.nc.sc <- integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
              integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)




# after bivariate transformation





bivariate.test3 <- function(x,y){
  exp(
      log(pgnorm(q = 1*(1-y > 1) + (1 - y)*(1 - y <= 1), mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -y*(y<=0),                     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)))
}

test3 <- function(x){
  exp(
    log(pgnorm(q = 1*(1-x > 1) + (1 - x)*(1 - x <= 1), mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
          pgnorm(q = -x*(x<=0),                     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)))
}
integral(test3, xmin = -1, xmax = 1) -
(integral(test3, xmin = -1, xmax = 0) +
integral(test3, xmin = 0,  xmax = 1))


integrate_debug(bivariate.test3, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = function(x) 1 - x)
integrate_debug(bivariate.test3, xmin = 0, xmax = 1, ymin = function(x) -x, ymax = 0) + 
  integrate_debug(bivariate.test3, xmin = 0, xmax = 1, ymin = 0,              ymax = function(x) 1 - x)

