rm(list=ls())

n1   <- 10
x1   <- 5
n2   <- 10
x2   <- 2

lik.1 <- function(eta1, eta2){
  dbinom(x1, size = n1, prob = eta1)*dbinom(x2, size = n2, prob = eta2)
}

integral2(lik.1, 0, 1, 0, 1)[[1]]

lik.2 <- function(theta, eta2){
  dbinom(x1,   size = n1, prob = theta + eta2)*dbinom(x2, size = n2, prob = eta2)
}

integral2(lik.2, -1, 0, function(theta) -theta, 1)[[1]] +
  integral2(lik.2, 0, 1, 0, function(theta) 1-theta)[[1]]


