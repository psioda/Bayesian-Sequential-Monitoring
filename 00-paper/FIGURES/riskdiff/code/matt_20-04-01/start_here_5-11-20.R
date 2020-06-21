evans2d <- function(n1, x1, m1, a1, b1, 
                    n2, x2, m2, a2, b2){
  
  post.nc       <- matrix(NA, nrow = n1 + 1, ncol = n2 + 1)
  
  
  for (i in 1:(n1 + 1)){
    print(i)
    for (j in 1:(n2 + 1)){
      post.kernel.1   <- function(x, y){
        exp(
          dbinom(i-1, n1, x + y, log = TRUE) +
          dbinom(j-1, n2, y, log = TRUE) + 
            dgnorm(x, m1, a1, b1, log = TRUE) - log(pgnorm(q = 1,   m1, a1, b1) - pgnorm(q = -1, m1, a1, b1)) +
            dgnorm(y, m2, a2, b2, log = TRUE) - log(pgnorm(q = 1-x, m2, a2, b2) - pgnorm(q = 0,  m2, a2, b2))
        )
      }
      post.kernel.2   <- function(x, y){
        exp(
          dbinom(i-1, n1, x + y, log = TRUE) +
          dbinom(j-1, n2, y, log = TRUE) +
            dgnorm(x, m1, a1, b1, log = TRUE) - log(pgnorm(q = 1, m1, a1, b1) - pgnorm(q = -1, m1, a1, b1)) +
            dgnorm(y, m2, a2, b2, log = TRUE) - log(pgnorm(q = 1, m2, a2, b2) - pgnorm(q = -x, m2, a2, b2))
        )
      }
      post.nc[i, j] <- integrate_debug(post.kernel.1, xmin = 0,  xmax = 1, ymin = 0,              ymax = function(x) 1 - x) +
                       integrate_debug(post.kernel.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
    }
  }
  print(sum(post.nc))
  sum(post.nc[post.nc <= post.nc[x1 + 1, x2 + 1]])
}


evans2d(y1.IP + y0.IP, y1.IP, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0,
        y1.PC + y0.PC, y1.PC, mu,         skpt.alpha0,    skpt.beta0)

, # always test at modal outcome
        9, 3, 1/3, 0.15, 2)

evans2d(9, 6, 1/3, 1, 2, # test unlikely point
        9, 3, 2/3, 1, 2)


enth.post.sc.1 <- function(x, y){
  dbinom(y1.IP, y1.IP + y0.IP, x + y)*
    dbinom(y1.PC, y1.PC + y0.PC, y)*
    dgnorm(x,         delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu,enth.alpha0, enth.beta0)/
    (pgnorm(q = 1 - x, mu, enth.alpha0, enth.beta0) -
       pgnorm(q = 0,   mu, enth.alpha0, enth.beta0))
}
enth.post.sc.2 <- function(x, y){
  dbinom(y1.IP, y1.IP + y0.IP, x + y)*
    dbinom(y1.PC, y1.PC + y0.PC, y)*
    dgnorm(x,         delta.enth, enth.rd.alpha0, enth.rd.beta0)/
    (pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) -
       pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0))*
    dgnorm(y,          mu,enth.alpha0, enth.beta0)/
    (pgnorm(q = 1,     mu, enth.alpha0, enth.beta0) -
       pgnorm(q = -x,  mu, enth.alpha0, enth.beta0))
}

# scaled normalizing constant for posterior density
# No longer scaled on 5/7/2020
skpt.nc.sc <- integrate_debug(skpt.post.sc.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
  integrate_debug(skpt.post.sc.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)