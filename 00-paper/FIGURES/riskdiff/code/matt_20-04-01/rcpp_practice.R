library(Rcpp)
library(bench)
library(RcppEigen)
library(RcppNumerical)
library(mvtnorm)

sourceCpp("2d_integral.cpp")
integrate_test2()

trueval = pmvnorm(c(-1, -1), c(1, 1), sigma = matrix(c(1, 0.5, 0.5, 1), 2))
as.numeric(trueval) - integrate_test2()$approximate








integrate_test()

sourceCpp("test1.cpp")


# 25.2 Getting started with C++


cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
            return sum;
            }')
# add works like a regular R function
add
#> function (x, y, z) 
#> .Call(<pointer: 0x7f08dc8450d0>, x, y, z)
add(1, 2, 3)

# 25.2.1 No inputs, scalar output
one <- function() 1L
one()

cppFunction('int one() {
  return 1;
}')

# 25.2.2 Scalar input, scalar output
cppFunction('int signC(int x) {
            if (x > 0) {
            return 1;
            } else if (x == 0) {
            return 0;
            } else {
            return -1;
            }
            }')
signC(-35)
signC(pi)

# 25.2.3 Vector input, scalar output
dat <- runif(4)

sumR <- function(x) {
  total <- 0
  for (i in seq_along(x)) {
    total <- total + x[i]
  }
  total
}
sumR(dat)

cppFunction('double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}')
sumC(dat)

x <- runif(1e3)
bench::mark(
  sum(x),
  sumC(x),
  sumR(x)
)[1:6]
