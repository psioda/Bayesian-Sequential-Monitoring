rm(list = ls())

#################################################################################################
## INITIAL CONDITIONS ###########################################################################
#################################################################################################

# lower target
theta_L<-0.15

# upper target
theta_H<-0.25

# tail probabilities for priors (low, high)
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

#################################################################################################
## INITIAL CONDITIONS ###########################################################################
#################################################################################################

## Prior model probabilities
prior_L<-1/2 
prior_H<-1/2 


#################################################################################################
## FREQUENTIST SAMPLE SIZE ######################################################################
#################################################################################################

freq_ss<-1000

#################################################################################################
## SIMULATIONS ##################################################################################
#################################################################################################

check_vector<-c(1,4,25,100,1000) # frequency of sequential monitoring
pi_range<-seq(theta_L,theta_H,by=0.01) # value of true response proportion
outer_trial_result_binary<-matrix(nrow=length(check_vector),ncol=length(pi_range)) # results matrix
outer_trial_result<-matrix(nrow=length(check_vector),ncol=length(pi_range)) # results matrix
outer_trial_result_ss<-matrix(nrow=length(check_vector),ncol=length(pi_range)) # results matrix

for (k in 1:length(check_vector)){

## OUTERMOST LOOP OVER TRUE PARAMETER VALUE ##
#outer_trial_result<-vector(length=length(pi_range))
#outer_trial_result_binary<-vector(length=length(pi_range))
#outer_freq_trial_result<-vector(length=length(pi_range))
#outer_trial_result_ss<-vector(length=length(pi_range)) # expected sample size
#outer_p_hat<-vector(length=length(pi_range))

for (j in 1:length(pi_range)){
  
  ## OUTER LOOP OVER TRIALS
  reps<-1000
  trial_result<-vector(length=reps)
  trial_result_binary<-vector(length=reps)
  trial_result_ss<-vector(length=reps)
  #freq_trial_result<-vector(length=reps)
  #inner_p_hat<-vector(length=reps)
  
  for (i in 1:reps){
    
    ## INNER LOOP OVER SAMPLE 
    futility<-0 # initialize
    efficacy<-0 # initialize

    n<-1 # initializing sample size
    y1_before<-0 # initializing number of successes
    
    while (n<=freq_ss){
      
      # generate trial result
      y<-rbinom(1,size=1,prob=pi_range[j])
      # new number of successes
      y1<-sum(y1_before)+y
      # new number of failures
      y0<-n-y1
      
      # futility (even optimst would give up)
      futility<-pbeta(theta_L,alpha_H+y1,beta_H+y0,lower.tail=TRUE)
      # efficacy (even pessimist would accept)
      efficacy<-pbeta(theta_H,alpha_L+y1,beta_L+y0,lower.tail=FALSE)

      ## for next loop
      y1_before<-y1
      n<-n+1
      
      if (n%%check_vector[k]==0 & (futility>.5 | efficacy>.5)){
        break
      }
    }
    
    ## Posterior model probabilities, slide 133/1005 of 779 notes
    # posterior probability of data given models
    post_L_D<-(beta(alpha_L,beta_L))^(-1)*beta(alpha_L+y1,beta_L+y0)
    post_H_D<-(beta(alpha_H,beta_H))^(-1)*beta(alpha_H+y1,beta_H+y0)
    # posterior probability of models given data
    post_L<-post_L_D*prior_L/(post_L_D*prior_L+post_H_D*prior_H)
    post_H<-post_H_D*prior_H/(post_L_D*prior_L+post_H_D*prior_H)
    ## compute probability of event A given the models, slide 7/1005 of 779 notes
    post_L_eventA<-pbeta(theta_L,alpha_L+y1,beta_L+y0,lower.tail=FALSE)
    post_H_eventA<-pbeta(theta_L,alpha_H+y1,beta_H+y0,lower.tail=FALSE)
    # compute probability of event A
    trial_result[i]<-sum(post_L*post_L_eventA,post_H*post_H_eventA)
    trial_result_binary[i]=(trial_result[i]>=0.95) # rejecting null hypothesis
    trial_result_ss[i]<-n-1
    #inner_p_hat[i]<-(y1/n-pi_range[j])
    #freq_trial_result[i]<-((sqrt(n)*(pi_range[j]-p_hat-delta_0)/sqrt(p_hat*(1-p_hat))<qnorm(alpha)) &
    #    (sqrt(n)*(pi_range[j]-p_hat+delta_0)/sqrt(p_hat*(1-p_hat))>qnorm(1-alpha)))
  }
  outer_trial_result_binary[k,j]<-mean(trial_result_binary) # sent to final results matrix
  outer_trial_result[k,j]<-mean(trial_result)
  #outer_freq_trial_result[j]<-mean(freq_trial_result)
  outer_trial_result_ss[k,j]<-mean(trial_result_ss)
  #outer_p_hat[j]<-mean(inner_p_hat)
}
}

outer_trial_result_binary # probability of rejecting null hypothesis (alpha at theta_L, 1-beta at theta_H)
outer_trial_result
outer_trial_result_ss


plot(pi_range,outer_trial_result_binary[1,],type='l')
lines(pi_range,outer_trial_result_binary[2,],type='l')
lines(pi_range,outer_trial_result_binary[3,],type='l')
lines(pi_range,outer_trial_result_binary[4,],type='l')


,dbeta(x,alpha_L,beta_L),type="l",col="red",xlab="x",ylab="f(x)",
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















outer_freq_trial_result

outer_p_hat


plot(p_range,power_plot,xlab="pi",ylab="power",type='l')
lines(pi_range,outer_trial_result_binary,xlab="pi",ylab="power",type='l',lty=5)
title(main="Bayesian power (dotted) vs Frequentist power (solid)")

plot(pi_range,outer_trial_result_ss,ylim=c(150,280),type='l',xlab='True pi',ylab='Expected sample size')

power_plot[1]
outer_trial_result_binary[1]




lines(pi_range,ss_100)
lines(pi_range,ss_20)
abline(h=271)
text(locator(), labels = c('No interm','Every 100','Every 20','Every subject'))

outer_p_hat
plot(pi_range,outer_p_hat,xlab='True p',ylab='bias')
#################################################################################################
#################################################################################################
#################################################################################################
