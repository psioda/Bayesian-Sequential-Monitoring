#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################
# Step 1: Create grid for possible values of phi
phi.range<-seq(0,100,by=0.001)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.skpt
quantiles.skpt<-qbeta(tail.skpt,(p.skpt)*phi.range,(1-(p.skpt))*phi.range,
                      lower.tail=FALSE)
# lower tail probability equal to tail.enth
quantiles.enth<-qbeta(tail.enth,(p.enth)*phi.range,(1-(p.enth))*phi.range,
                      lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi.range[which.min(abs(p.enth-quantiles.skpt))] # fixed 5/13/19
phi_H<-phi.range[which.min(abs(p.skpt-quantiles.enth))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.skpt<-(p.skpt)*phi_L
beta.skpt<-(1-(p.skpt))*phi_L
alpha.enth<-(p.enth)*phi_H
beta.enth<-(1-(p.enth))*phi_H

# Spike and slab prior information
mu0.skpt<-0.2
## flat
# sigma0.skpt<-0.2128
# lambda0.skpt<-4
sigma0.skpt<-0.1529
lambda0.skpt<-1.7

mu0.enth<-0.4
## flat
#sigma0.enth<-0.2149
#lambda0.enth<-4
## regular
sigma0.enth<-0.1906
lambda0.enth<-2.5