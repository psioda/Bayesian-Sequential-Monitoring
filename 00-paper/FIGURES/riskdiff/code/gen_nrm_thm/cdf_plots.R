rm(list = ls())
library(gnorm)

x <- seq(-10, 10, length = 1000)
beta <- seq(2, 10, length = 20)
alpha <- seq(0.1, 100, length = 20)

y <- pgnorm(x, m = 0, alpha = 1, beta = 2)
plot(x, y, type = 'l')

for (i in 1:length(beta)){
  for (j in 1:length(alpha)){
  y <- pgnorm(x, m = 0, alpha = alpha[j], beta = beta[i])
  lines(x, y, type = 'l')
  }
}