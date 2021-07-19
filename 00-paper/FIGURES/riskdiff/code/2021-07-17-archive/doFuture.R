library(doFuture)
registerDoFuture()
plan(multiprocess)

Nlength <- 1000
Nvector <- 3

res <- foreach(i = 1:Nvector, .combine = cbind) %dopar% {
  1:Nlength / pi * rnorm(1)
}
res <- data.frame(res)
colnames(res) <- paste0("vector", 1:Nvector)

dim(res)