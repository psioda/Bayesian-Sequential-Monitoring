# Name: Evan Kwiatkowski
# Date: April 1, 2020
# Purpose: Explore adaptive weight for skeptical monitoring prior 

library(gnorm)

## parameters for joint distribution of theta (risk difference) and psi (placebo response)
## pi(theta,psi)=pi(theta)*pi(psi|theta)
## pi(theta),pi(psi|theta) ~ generalizedNormal()
# enthusiastic parameters
theta.e <- c(0.12, 0.08658609, 2)        # mu, alpha, beta
psi.e   <- c(0.39, 0.2190863,  1.944127) # mu, alpha, beta
# skeptical parameters
theta.s <- c(0,    0.05776004, 1.282048) # mu, alpha, beta
psi.s   <- c(0.39, 0.2158015,  1.920246) # mu, alpha, beta

# source priors
# note prior.1 is for nonnegative risk difference
# note prior.2 is for negative risk difference
source("priors.R")

# specify the observed risk difference and placebo response probability
rd.mle <- seq(0, 0.12, by = 0.01)
pc.mle <- 0.39

# evaluate and print the weight for skeptical monitoring prior 
for (i in 1:length(rd.mle)){
# likelihoods
if (rd.mle[i] >= 0){
  skpt.lik <- skpt.prior.1(rd.mle[i], pc.mle)
  enth.lik <- enth.prior.1(rd.mle[i], pc.mle)
} else {
  skpt.lik <- skpt.prior.2(rd.mle[i], pc.mle)
  enth.lik <- enth.prior.2(rd.mle[i], pc.mle)
}

skpt.mix.wt <- skpt.lik/sum(skpt.lik, enth.lik)

print(paste0("risk difference: ", rd.mle[i], " skeptical mixture weight:", skpt.mix.wt))

}
