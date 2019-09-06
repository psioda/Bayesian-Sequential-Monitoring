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

# ## Design parameters, default case
# alpha.skpt.1<-alpha.skpt
# beta.skpt.1<-beta.skpt
# alpha.skpt.2<-alpha.skpt
# beta.skpt.2<-beta.skpt
# mix.1<-0.5
# alpha.enth.1<-alpha.enth
# beta.enth.1<-beta.enth
# alpha.enth.2<-alpha.enth
# beta.enth.2<-beta.enth
# mix.2<-0.5
# ## Design parameters, default skeptical, spike/slab enthuastic
# alpha.skpt.1<-alpha.skpt
# beta.skpt.1<-beta.skpt
# alpha.skpt.2<-alpha.skpt
# beta.skpt.2<-beta.skpt
# mix.1<-0.5
# alpha.enth.1<-2.460314
# beta.enth.1<-3.603074
# alpha.enth.2<-96.979737
# beta.enth.2<-143.682052
# mix.2<-0.3426655
# 
# ## Design parameters, default enthuastic, spike/slab skeptical
# alpha.skpt.1<-0.969595
# beta.skpt.1<-2.186593
# alpha.skpt.2<-10.199189
# beta.skpt.2<-46.344410
# mix.1<-0.1421888
# alpha.enth.1<-alpha.enth
# beta.enth.1<-beta.enth
# alpha.enth.2<-alpha.enth
# beta.enth.2<-beta.enth
# mix.2<-0.5

# ## Design parameters, spike/slab case ##
# alpha.skpt.1<-0.969595  
# beta.skpt.1<-2.186593 
# alpha.skpt.2<-10.199189 
# beta.skpt.2<-46.344410 
# mix.1<-0.1421888
# alpha.enth.1<-2.460314 
# beta.enth.1<-3.603074 
# alpha.enth.2<-96.979737
# beta.enth.2<-143.682052
# mix.2<-0.3426655 