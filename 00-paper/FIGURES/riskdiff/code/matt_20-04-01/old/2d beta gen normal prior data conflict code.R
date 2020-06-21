rm(list=ls())
library(pracma)
library(gnorm)
evans2d <- function(n1, x1, m1, a1, b1, 
                    n2, x2, m2, a2, b2){
  
  post.nc       <- matrix(NA, nrow = n1 + 1, ncol = n2 + 1)
  
  
  for (i in 1:(n1+1)){
    for (j in 1:(n2+1)){
      post.kernel   <- function(x, y) 
        exp(
          dbinom(i-1, n1, x, log = TRUE) +
            dbinom(j-1, n2, y, log = TRUE) + 
            dgnorm(x, m1, a1, b1, log = TRUE) - log(pgnorm(q = 1, m1, a1, b1) - pgnorm(q = 0, m1, a1, b1)) +
            dgnorm(y, m2, a2, b2, log = TRUE) - log(pgnorm(q = 1, m2, a2, b2) - pgnorm(q = 0, m2, a2, b2))
        )
      post.nc[i, j] <- integral2(post.kernel, 0, 1, 0, 1, singular = TRUE)[[1]]
    }
  }
  sum(post.nc[post.nc <= post.nc[x1 + 1, x2 + 1]])
}


evans2d(9, 3, 1/3, 1, 2, # always test at modal outcome
        9, 6, 2/3, 1, 2)

evans2d(9, 6, 1/3, 1, 2, # test unlikely point
        9, 3, 2/3, 1, 2)
