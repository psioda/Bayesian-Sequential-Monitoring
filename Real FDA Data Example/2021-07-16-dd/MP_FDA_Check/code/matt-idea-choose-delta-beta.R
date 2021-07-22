rm(list = ls())

theta0 <- 0.45
theta1 <- 0.55
theta2 = 0.65

phis <- 1
tails <- 1
while (tails >= 0.025){
  phis <- phis + 0.001
  tails <- pbeta(theta1,theta0*phis,(1-theta0)*phis, lower.tail = F)
}

phie <- 1
taile <- 1
while (taile >= 0.025){
  phie <- phie + 0.001
  taile <- pbeta(theta0,theta1*phie,(1-theta1)*phie, lower.tail = T)
}

y.tot <- 50
y.seq <- seq(0, y.tot)
BoxPE <- BoxPS <- rep(NA, length(y.seq))
for (i in 1:length(y.seq)){
  BoxPS[i] <- sum(dbbinom(y.seq, y.tot,theta0*phis,(1-theta0)*phis) * 
                    (dbbinom(y.seq, y.tot,theta0*phis,(1-theta0)*phis) <= dbbinom(i, y.tot,theta0*phis,(1-theta0)*phis)))
  
  BoxPE[i] <- sum(dbbinom(y.seq, y.tot,theta1*phie,(1-theta1)*phie) * 
                    (dbbinom(y.seq, y.tot,theta1*phie,(1-theta1)*phie) <= dbbinom(i, y.tot,theta1*phie,(1-theta1)*phie)))
}


plot(y.seq, BoxPS, type = 'l')
lines(y.seq, BoxPE)

delta <- seq(0, 0.25, by = 0.05)
betas <- seq(1, 5, by = 0.5)
for (i in 1:length(delta)){
  for (j in 1:length(betas)){
    lines(y.seq, (1 - delta[i])*pbeta(BoxPE, 1, betas[j]))
  }
}



x.seq <- seq(0, 1, length = 1E3 + 1)
x.tra <- pbeta(x.seq, 1, 5)
plot(x.seq, x.tra, type = 'l')
abline(v = 0.5)
x.tra[x.seq == 0.5]
x.tra[x.seq == 0.05]
x.tra[x.seq == 0.025]


d.list <- seq(0.5, 0.9, by = 0.1)
b.list <- seq(1, 5, length = 4E4 + 1)
final <- rep(NA, length(d.list))
for (i in 1:length(d.list)){
  ans <- 0.5
  j <- 1
  while(ans <= d.list[i]){
    j <- j + 1
    ans <- pbeta(0.5, 1, b.list[j])
  }
  final[i] <- b.list[j]
}




