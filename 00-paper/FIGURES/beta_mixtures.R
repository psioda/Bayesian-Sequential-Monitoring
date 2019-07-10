rm(list = ls())

#################################################################################################
## INITIAL CONDITIONS ###########################################################################
#################################################################################################

# lower target
theta_L<-0.15

# upper target
theta_H<-0.45

# tail probabilities for priors (low, middle, high)
alpha_L<-0.05
alpha_H<-0.05


#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################

## Parameterize prior models
# Step 1: Create grid for possible values of phi
# E.g. beta parameterization for middle prior: beta(pi_0*phi,(1-pi_0)*phi)
phi_seq<-seq(0,100,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to alpha_L/2
test_L<-qbeta(alpha_L/2,(theta_L)*phi_seq,(1-(theta_L))*phi_seq,lower.tail=FALSE)
# lower tail probability equal to alpha_H/2
test_H<-qbeta(alpha_H/2,(theta_H)*phi_seq,(1-(theta_H))*phi_seq,lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi_seq[which.min(abs(theta_H-test_L))] # fixed 5/13/19
phi_H<-phi_seq[which.min(abs(theta_L-test_H))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha_L<-(theta_L)*phi_L
beta_L<-(1-(theta_L))*phi_L

alpha_H<-(theta_H)*phi_H
beta_H<-(1-(theta_H))*phi_H

# Plot results
x<-seq(0,1,by=0.01)
#low (skeptical)
plot(x,dbeta(x,alpha_L,beta_L),type="l",col="red",xlab="x",ylab="f(x)",
     ylim=c(0,max(dbeta(x,alpha_L,beta_L),dbeta(x,alpha_H,beta_H))))

#high (enthuastic)
lines(x,dbeta(x,alpha_H,beta_H),type="l",col="blue")
points(x,1/4*dbeta(x,alpha_L,beta_L)+3/4*dbeta(x,alpha_H,beta_H),type='l',lty=5,lwd=2)
points(x,1/2*dbeta(x,alpha_L,beta_L)+1/2*dbeta(x,alpha_H,beta_H),type='l',lty=6,lwd=2)
points(x,3/4*dbeta(x,alpha_L,beta_L)+1/4*dbeta(x,alpha_H,beta_H),type='l',lty=7,lwd=2)
title(main="Beta priors")
title(main="Beta priors")
abline(v=theta_L)
abline(v=theta_H)

#abline(v=(pi_0-delta_0))
#abline(v=pi_0)
#abline(v=(pi_0+delta_0))
