skpt.post.1.b <- function(x, y){ # for x > 0 (theta > 0)
  dbinom(y1.IP, y1.IP + y0.IP, x + y)*
  dbinom(y1.PC, y1.PC + y0.PC, y)*
    dgnorm(x,         delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
    (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
       pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
    dgnorm(y,          mu,skpt.alpha0, skpt.beta0)/
    (pgnorm(q = 1 - x, mu, skpt.alpha0, skpt.beta0) -
       pgnorm(q = 0,   mu, skpt.alpha0, skpt.beta0))
}
skpt.post.2.b <- function(x, y){ # for x < 0 (theta < 0)
  dbinom(y1.IP, y1.IP + y0.IP, x + y)*
  dbinom(y1.PC, y1.PC + y0.PC, y)*
    dgnorm(x,         delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)/
    (pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) -
       pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0))*
    dgnorm(y,          mu,skpt.alpha0, skpt.beta0)/
    (pgnorm(q = 1,     mu, skpt.alpha0, skpt.beta0) -
       pgnorm(q = -x,  mu, skpt.alpha0, skpt.beta0))
}
enth.post.1.b <- function(x, y){ # for x > 0 (theta > 0)
  dbinom(y1.IP, y1.IP + y0.IP, x + y)*
  dbinom(y1.PC, y1.PC + y0.PC, y)*
    dgnorm(x,         delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu,enth.alpha0, enth.beta0)/
    (pgnorm(q = 1 - x, mu, enth.alpha0, enth.beta0) -
       pgnorm(q = 0,   mu, enth.alpha0, enth.beta0))
}
enth.post.2.b <- function(x, y){ # for x < 0 (theta < 0)
  dbinom(y1.IP, y1.IP + y0.IP, x + y)*
  dbinom(y1.PC, y1.PC + y0.PC, y)*
    dgnorm(x,         delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu,enth.alpha0, enth.beta0)/
    (pgnorm(q = 1,     mu, enth.alpha0, enth.beta0) -
       pgnorm(q = -x,  mu, enth.alpha0, enth.beta0))
}

skpt.nc.sc.1.b <- integrate_debug(skpt.post.1.b, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x)
skpt.nc.sc.2.b <- integrate_debug(skpt.post.2.b, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
enth.nc.sc.1.b <- integrate_debug(enth.post.1.b, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x)
enth.nc.sc.2.b <- integrate_debug(enth.post.2.b, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)

skpt.nc.sc.b <- skpt.nc.sc.1.b + skpt.nc.sc.2.b
enth.nc.sc.b <- enth.nc.sc.1.b + enth.nc.sc.2.b
fut.skpt.wt.b <- fut.mix.prob*skpt.nc.sc.b/(fut.mix.prob*skpt.nc.sc.b + (1 - fut.mix.prob)*enth.nc.sc.b)
fut.skpt.wt.b
# 
# # un-normalized posterior density
# skpt.post.1 <- function(x, y){
#   exp(
#     y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
#       y1.PC*log(y)     + y0.PC*log(1 - y) -
#       (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
#       (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
#       log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
#             pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
#   )
# }
# skpt.post.2 <- function(x, y){
#   exp(
#     y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
#       y1.PC*log(y)     + y0.PC*log(1 - y) -
#       (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
#       (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
#       log(pgnorm(q = 1,  mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
#             pgnorm(q = -x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0))
#   )
# }
# enth.post.1 <- function(x, y){
#   exp(
#     y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
#       y1.PC*log(y)     + y0.PC*log(1 - y) -
#       (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
#       (abs(y - mu)/enth.alpha0)^enth.beta0 -
#       log(pgnorm(q = 1 - x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
#             pgnorm(q = 0,     mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
#   )
# }
# enth.post.2 <- function(x, y){
#   exp(
#     y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
#       y1.PC*log(y)     + y0.PC*log(1 - y) -
#       (abs(x - delta.enth)/enth.rd.alpha0)^enth.rd.beta0 -
#       (abs(y - mu)/enth.alpha0)^enth.beta0 -
#       log(pgnorm(q = 1,  mu = mu, alpha = enth.alpha0, beta  = enth.beta0) -
#             pgnorm(q = -x, mu = mu, alpha = enth.alpha0, beta  = enth.beta0))
#   )
# }
# 
# # scaled normalizing constant for posterior density
# skpt.nc.sc <- integrate_debug(skpt.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#   integrate_debug(skpt.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
# enth.nc.sc <- integrate_debug(enth.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
#   integrate_debug(enth.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
# 
# # posterior mixing weights
# # http://www.mas.ncl.ac.uk/~nmf16/teaching/mas3301/week11.pdf
# fut.skpt.wt <- fut.mix.prob*skpt.nc.sc/(fut.mix.prob*skpt.nc.sc + (1 - fut.mix.prob)*enth.nc.sc)
# 
# 
# function(p1, p2) dbinom(i-1, n1, p1)*
#                                         dbinom(j-1, n2, p2)*
#                                         dbeta(p1, a1, b1)*
#                                         dbeta(p2, a2, b2)
#       #post.nc <- integral2(post.kernel, 1E-1, 1, 1E-1, 1, singular = TRUE)[[1]]
#       post.nc[i, j] <- integral2(post.kernel, 0, 1, 0, 1, singular = TRUE)[[1]]
#     }
#   }
#   sum(post.nc[post.nc <= post.nc[x1 + 1, x2 + 1]])
# }
# evans2d_log <- function(n1, x1, a1, b1, 
#                         n2, x2, a2, b2){
#   
#   post.nc       <- matrix(NA, nrow = n1 + 1, ncol = n2 + 1)
#   
#   
#   for (i in 1:(n1+1)){
#     for (j in 1:(n2+1)){
#       #post.nc <- NA
#       post.kernel   <- function(p1, p2) 
#         dbinom(i-1, n1, p1, log = TRUE) +
#         dbinom(j-1, n2, p2, log = TRUE) +
#         dbeta(p1, a1, b1, log = TRUE) +
#         dbeta(p2, a2, b2, log = TRUE)
#       #post.nc <- integral2(post.kernel, 1E-1, 1, 1E-1, 1, singular = TRUE)[[1]]
#       post.nc[i, j] <- integral2(post.kernel, 0, 1, 0, 1, singular = TRUE)[[1]]
#     }
#   }
#   sum(post.nc[post.nc <= post.nc[x1 + 1, x2 + 1]])
# }
# 
# 
# # seems to be working
# evans2d(10,3,4,8,20,12,8,4)
# evans2d(10,8,4,8,20,4,8,4)
# # next use post.nc to verify a value...
# 
# 
# 
# 
# # goal: work up to showing that mixing weight (things from nc.blah blah) are the same
# # will need to see if R understands using Rbinom with (x+y) as response probability, then integrating
# 
# # is not the end of the world if it doesn't work
# 
# # or could skip this step for now ... 
# 
# 
# # big picture anything that shows the 2d prior predictive distribution and prior data conflict calculation are working
# 
# evans <- function(n, t0, a, b){
#   
#   n<-10
#   a<-4
#   b<-8
#   t0<-5
#   
#   prior.kernel <- function(p) p^(a-1)*(1-p)^(b-1)
#   prior.nc     <- integrate(prior.kernel, 0, 1)[[1]]
#   
#   
#   post.nc      <- NA
#   for (x in 1:(n+1)){ # has to start binomial looping at 0
#     post.kernel   <- function(p) choose(n,(x-1))*p^((x-1)+a-1)*(1-p)^(n-(x-1)+b-1)/prior.nc
#     post.nc[x]    <- integrate(post.kernel, 0, 1)[[1]]
#   }
#   
#   sum(post.nc[post.nc <= post.nc[t0 + 1]])
#   
# }
# 
# 
# # scaled (un-normalized) posterior density
# skpt.post.sc.1 <- function(x, y){
#   exp(
#     y1.IP*log(x + y) + y0.IP*log(1 - (x + y)) + 
#       y1.PC*log(y)     + y0.PC*log(1 - y) -
#       (abs(x - delta.skpt)/skpt.rd.alpha0)^skpt.rd.beta0 -
#       (abs(y - mu)/skpt.alpha0)^skpt.beta0 -
#       log(pgnorm(q = 1 - x, mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0) -
#             pgnorm(q = 0,     mu = mu, alpha = skpt.alpha0, beta  = skpt.beta0)) -
#       sc
#   )
# }
# 
# 
# test <- function(p1,p2){
#   dbinom(x, n, p1+p2)
# }
# n<-10
# x<-5
# library("pracma")
# integrate_debug(test,0,1,0,0.5)
# 
