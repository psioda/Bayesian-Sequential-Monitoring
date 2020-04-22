rm(list = ls())
library("rmutil")

ss      <- 15
alpha_e <- 8
beta_e  <- 4
alpha_s <- 4
beta_s  <- 8
y_obs   <- 10
  
y   <- seq(0, ss)

d_e   <- dbetabinom(y = y,
                    size = ss,
                    m = alpha_e/(alpha_e+beta_e),
                    s = alpha_e+beta_e)
d_s   <- dbetabinom(y = y,
                    size = ss,
                    m = alpha_s/(alpha_s+beta_s),
                    s = alpha_s+beta_s)

plot(y, d_e)
plot(y, d_s)

d_e_obs <- dbetabinom(y = y_obs,
                      size = ss,
                      m = alpha_e/(alpha_e+beta_e),
                      s = alpha_e+beta_e)
d_s_obs <- dbetabinom(y = y_obs,
                      size = ss,
                      m = alpha_s/(alpha_s+beta_s),
                      s = alpha_s+beta_s)

sum(d_e*(d_e <= d_e_obs))
sum(d_s*(d_s <= d_s_obs))

pred_e <- NA
pred_s <- NA

for (y in 0:ss){
  idx <- y + 1
  pred_e[idx] <- (choose(ss,y)/beta(alpha_e,beta_e))*integrate(function(p) (p^(y+alpha_e-1))*((1-p)^(ss-y+beta_e-1)), 0, 1)$value
  pred_s[idx] <- (choose(ss,y)/beta(alpha_s,beta_s))*integrate(function(p) (p^(y+alpha_s-1))*((1-p)^(ss-y+beta_s-1)), 0, 1)$value
}
pred_e_obs <- (choose(ss,y_obs)/beta(alpha_e,beta_e))*integrate(function(p) (p^(y_obs+alpha_e-1))*((1-p)^(ss-y_obs+beta_e-1)), 0, 1)$value
pred_s_obs <- (choose(ss,y_obs)/beta(alpha_s,beta_s))*integrate(function(p) (p^(y_obs+alpha_s-1))*((1-p)^(ss-y_obs+beta_s-1)), 0, 1)$value

sum(pred_e)
sum(pred_s)

sum(pred_e*(pred_e <= pred_e_obs))
sum(pred_s*(pred_s <= pred_s_obs))
