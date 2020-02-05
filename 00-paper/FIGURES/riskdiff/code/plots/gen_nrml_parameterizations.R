alpha.seq <- seq(1E-2, 2, length = 1000)
beta.seq  <- seq(1E-2, 2, length = 1000)

result1 <- matrix(NA, nrow = length(alpha.seq), ncol = length(beta.seq))
result2 <- matrix(NA, nrow = length(alpha.seq), ncol = length(beta.seq))
  
q.outer <- 0.025
q.inner <- 0.125
  
theta0 <- 0
theta1 <- 1

for (i in 1:length(alpha.seq)){
  if (i%%10 == 0){print(paste0("Simulation ", i))}
  for (j in 1:length(beta.seq)){
    alpha <- alpha.seq[i]
    beta  <- beta.seq[j]
    result1[i,j] <- qgnorm(q.outer,
                           mu = theta1,
                           alpha = alpha,
                           beta = beta)
    result2[i,j] <-qgnorm(q.inner,
                          mu = theta1,
                          alpha = alpha,
                          beta = beta)
  }
}

result3 <- sqrt((result1 - theta0)^2) + sqrt((result2 - (theta0+theta1)/2)^2)

index <- which(result3  ==  min(result3), arr.ind = TRUE)
i <- index[1]
j <- index[2]

alpha <- alpha.seq[i]
beta  <- beta.seq[j]

# check answer
pgnorm(theta0, mu = theta1, alpha = alpha, beta = beta)
pgnorm((theta0+theta1)/2, mu = theta1, alpha = alpha, beta = beta)
alpha
beta

# q.outer <- 0.025
# q.inner <- 0.075
#   > pgnorm(theta0, mu = theta1, alpha = alpha, beta = beta)
# [1] 0.02503618
# > pgnorm((theta0+theta1)/2, mu = theta1, alpha = alpha, beta = beta)
# [1] 0.07498953
# > alpha
# [1] 0.03893894
# > beta
# [1] 0.4856757

# q.outer <- 0.025
# q.inner <- 0.175
# > # check answer
#   > pgnorm(theta0, mu = theta1, alpha = alpha, beta = beta)
# [1] 0.02499561
# > pgnorm((theta0+theta1)/2, mu = theta1, alpha = alpha, beta = beta)
# [1] 0.1750171
# > alpha
# [1] 0.7842843
# > beta
# [1] 2.306807

x <- seq(-3,3,length=1000)
y1 <- dgnorm(x, mu = theta1, alpha = 0.4502302, beta = 1.217147)
y2 <- dgnorm(x, mu = theta1, alpha = 0.7842843, beta = 2.306807)

plot(x,y1,type='l')
lines(x, y2)

plot(c(rd.grid.lower, rd.grid.upper),
     c(skpt.lower, skpt.upper),
     type = 'l',
     xlim = c(-0.4, 0.4),
     ylim = c(0, 9.25),
     xlab = "",
     ylab = "",
     main = "",
     xaxt = "n",
     yaxt = "n")

title(ylab = "Density Value", 
      line = 1)

title(xlab = "Response Probability",
      line = 2)

axis(1,
     at = c(delta.skpt,delta.enth),
     labels = c(as.expression(bquote(delta[S])),
                as.expression(bquote(delta[E]))))

polygon(c(rd.grid.lower, rd.grid.upper[rd.grid.upper <= delta.enth], delta.enth),
        c(skpt.lower, skpt.upper[rd.grid.upper <= delta.enth], 0),
        col = "lightgrey")

polygon(c(rd.grid.upper[rd.grid.upper >= delta.enth], delta.enth),
        c(skpt.upper[rd.grid.upper >= delta.enth], 0),col = "black")

segments(x0 = delta.skpt, 
         y0 = 0, 
         y1 = skpt.upper[which(rd.grid.upper == delta.skpt)])
segments(x0 = delta.enth, 
         y0 = 0, 
         y1 = skpt.upper[which(rd.grid.upper == delta.enth)])

legend('topright',
       legend= c(as.expression(bquote(P(theta[1]-theta[0]>delta[E])==.(1-sig.eff)))))
