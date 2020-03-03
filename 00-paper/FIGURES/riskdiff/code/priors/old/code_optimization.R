rm(list = ls())

mu <- 0

p1 <- 0.95
q1 <- 3
p2 <- 0.90
q2 <- 2

f <- function(x){
 (pgnorm(q1, mu = mu, alpha = x[1], beta = x[2]) - p1)^2 +
 (pgnorm(q2, mu = mu, alpha = x[1], beta = x[2]) - p2)^2
}

optim(c(1,1), f)
x <- optim(c(1,1), f)$par

# check
pgnorm(q1, mu = mu, alpha = x[1], beta = x[2])
pgnorm(q2, mu = mu, alpha = x[1], beta = x[2])

# Derivative-free optimization methods of Nelder and Mead (1965)