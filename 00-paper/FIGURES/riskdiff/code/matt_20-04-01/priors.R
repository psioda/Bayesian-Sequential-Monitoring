# normalized prior densities
enth.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  exp(-(abs(x - theta.e[1])/theta.e[2])^theta.e[3])/(2*theta.e[2]*gamma(1/theta.e[3])/theta.e[3])*
  exp(-(abs(y - psi.e[1])/psi.e[2])^psi.e[3])/(2*psi.e[2]*gamma(1/psi.e[3])/psi.e[3])/
    (pgnorm(q = 1,     mu = theta.e[1], alpha = theta.e[2], beta = theta.e[3]) -
     pgnorm(q = -1,    mu = theta.e[1], alpha = theta.e[2], beta = theta.e[3]))/
    (pgnorm(q = 1 - x, mu = psi.e[1],   alpha = psi.e[2],   beta  = psi.e[3]) -
     pgnorm(q = 0,     mu = psi.e[1],   alpha = psi.e[2],   beta  = psi.e[3]))
}
enth.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  exp(-(abs(x - theta.e[1])/theta.e[2])^theta.e[3])/(2*theta.e[2]*gamma(1/theta.e[3])/theta.e[3])*
  exp(-(abs(y - psi.e[1])/psi.e[2])^psi.e[3])/(2*psi.e[2]*gamma(1/psi.e[3])/psi.e[3])/
    (pgnorm(q = 1,  mu = theta.e[1], alpha = theta.e[2], beta = theta.e[3]) -
     pgnorm(q = -1, mu = theta.e[1], alpha = theta.e[2], beta = theta.e[3]))/
    (pgnorm(q = 1,  mu = psi.e[1],   alpha = psi.e[2],   beta = psi.e[3]) -
     pgnorm(q = -x, mu = psi.e[1],   alpha = psi.e[2],   beta = psi.e[3]))
}

skpt.prior.1 <- function(x, y){ # for x > 0 (theta > 0)
  exp(-(abs(x - theta.s[1])/theta.s[2])^theta.s[3])/(2*theta.s[2]*gamma(1/theta.s[3])/theta.s[3])*
  exp(-(abs(y - psi.s[1])/psi.s[2])^psi.s[3])/(2*psi.s[2]*gamma(1/psi.s[3])/psi.s[3])/
    (pgnorm(q = 1,     mu = theta.s[1], alpha = theta.s[2], beta = theta.s[3]) -
     pgnorm(q = -1,    mu = theta.s[1], alpha = theta.s[2], beta = theta.s[3]))/
    (pgnorm(q = 1 - x, mu = psi.s[1],   alpha = psi.s[2],   beta  = psi.s[3]) -
     pgnorm(q = 0,     mu = psi.s[1],   alpha = psi.s[2],   beta  = psi.s[3]))
}
skpt.prior.2 <- function(x, y){ # for x < 0 (theta < 0)
  exp(-(abs(x - theta.s[1])/theta.s[2])^theta.s[3])/(2*theta.s[2]*gamma(1/theta.s[3])/theta.s[3])*
  exp(-(abs(y - psi.s[1])/psi.s[2])^psi.s[3])/(2*psi.s[2]*gamma(1/psi.s[3])/psi.s[3])/
    (pgnorm(q = 1,  mu = theta.s[1], alpha = theta.s[2], beta = theta.s[3]) -
     pgnorm(q = -1, mu = theta.s[1], alpha = theta.s[2], beta = theta.s[3]))/
    (pgnorm(q = 1,  mu = psi.s[1],   alpha = psi.s[2],   beta  = psi.s[3]) -
     pgnorm(q = -x, mu = psi.s[1],   alpha = psi.s[2],   beta  = psi.s[3]))
}