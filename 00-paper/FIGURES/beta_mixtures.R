rm(list = ls())

# mean of skeptical prior
theta_L<-0.20
# mean of enthuastic prior
theta_H<-0.40
# futility theta
theta_futility<-0.30
# value of true response proportion
pi_range<-seq(theta_L-0.05,theta_H+0.05,by=0.05)

# tail probabilities for priors (low, high)
alpha_L<-0.045
alpha_H<-0.05
# Prior model probabilities
prior_L<-1/2 
prior_H<-1/2 
# maximum sample sizes
max_ss<-76
# frequency of sequential monitoring
check_vector<-c(2) 
# significant trial result threshold
significant_threshold<-0.9
futility_threshold<-0.85
efficacy_threshold<-0.95
# things to do with futility & efficacy


#example1<-function(){

#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################

# Step 1: Create grid for possible values of phi
phi_seq<-seq(0,100,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to alpha_L
test_L<-qbeta(alpha_L,(theta_L)*phi_seq,(1-(theta_L))*phi_seq,lower.tail=FALSE)
# lower tail probability equal to alpha_H
test_H<-qbeta(alpha_H,(theta_H)*phi_seq,(1-(theta_H))*phi_seq,lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi_seq[which.min(abs(theta_H-test_L))] # fixed 5/13/19
phi_H<-phi_seq[which.min(abs(theta_L-test_H))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha_L<-(theta_L)*phi_L
beta_L<-(1-(theta_L))*phi_L
alpha_H<-(theta_H)*phi_H
beta_H<-(1-(theta_H))*phi_H

# Step 5: Plot Skeptical and enthuastic priors separately
## 7/19/19 match slide 14 from FDA presentation, shade in gray scale, set 350 DPI
par(ask=TRUE)
par(mfrow = c(1, 2)) 
x<-seq(0,1,by=0.01)
# Make 2 boxplots
#low (skeptical)
plot(x,dbeta(x,alpha_L,beta_L),type="l",col="red",xlab="Response Probability",ylab="Density Value",main="Skeptical Prior",
     #xaxt="n",
     ylim=c(0,max(dbeta(x,alpha_L,beta_L),dbeta(x,alpha_H,beta_H))))
#axis(1,at=theta_L,labels=expression(theta[0]))
#axis(1,at=seq(0,1,by=0.2))
abline(v=theta_H)
abline(v=theta_L)
#high (enthuastic)
plot(x,dbeta(x,alpha_H,beta_H),type="l",col="blue",xlab="Response Probability",ylab="Density Value",main="Enthuastic Prior",
     #xaxt="n",
     ylim=c(0,max(dbeta(x,alpha_L,beta_L),dbeta(x,alpha_H,beta_H))))
#axis(1,at=theta_H,labels=expression(theta[A]))
abline(v=theta_H)
abline(v=theta_L)
# Step 6: Plot inference priors (mixtures)
## 7/19/19 make 50:50 solid black line, 25:75 dark grey with different dash patterns,
# make skeptical/enthuastic light grey with smaller thinkness
# no legend, put omega on actual lines

par(mfrow = c(1, 1)) 
plot(x,dbeta(x,alpha_L,beta_L),type="l",col="red",xlab="Response Probability",ylab="Density Value",main="Inference Priors",
     xaxt="n",
     ylim=c(0,max(dbeta(x,alpha_L,beta_L),dbeta(x,alpha_H,beta_H))))
axis(1,at=theta_L,labels=expression(theta[0]))
lines(x,dbeta(x,alpha_H,beta_H),type="l",col="blue")
axis(1,at=theta_H,labels=expression(theta[A]))
axis(1,at=theta_L,labels=expression(theta[0]))
points(x,1/4*dbeta(x,alpha_L,beta_L)+3/4*dbeta(x,alpha_H,beta_H),type='l',lty=3)
points(x,1/2*dbeta(x,alpha_L,beta_L)+1/2*dbeta(x,alpha_H,beta_H),type='l',lty=2)
points(x,3/4*dbeta(x,alpha_L,beta_L)+1/4*dbeta(x,alpha_H,beta_H),type='l',lty=1)
legend(x=0.55,y=4.5,
       legend=c("Skeptical","75:25","50:50","25:75","Enthuastic"),
       col=c("red","black","black","black","blue"),lty=c(1,1,2,3,1),
       cex=0.8,text.font=3)

#################################################################################################
## SIMULATIONS ##################################################################################
#################################################################################################

## Step 1: Create outer loop based on frequency of interim analyses
# probability of proof-of-effect
outer_trial_result<-matrix(nrow=length(check_vector),ncol=length(pi_range)) 
# binary decision of proof-of-effect
outer_trial_result_binary<-matrix(nrow=length(check_vector),ncol=length(pi_range)) 
# overall sample size
outer_trial_result_ss<-matrix(nrow=length(check_vector),ncol=length(pi_range)) 
# sample mean
outer_p_hat<-matrix(nrow=length(check_vector),ncol=length(pi_range))
# stop early for efficacy
outer_efficacy<-matrix(nrow=length(check_vector),ncol=length(pi_range))
# stop early for futility
outer_futility<-matrix(nrow=length(check_vector),ncol=length(pi_range))
# inconclusive findings with full dataset
outer_inconclusive<-matrix(nrow=length(check_vector),ncol=length(pi_range))
# posterior mean
outer_posterior_mean<-matrix(nrow=length(check_vector),ncol=length(pi_range))
# coverage probability
outer_coverage<-matrix(nrow=length(check_vector),ncol=length(pi_range))


#outer_freq_trial_result<-vector(length=length(pi_range))

for (k in 1:length(check_vector)){

## Step 2: Create outer loop on range of true response values (for Type 1 error and power curves)
for (j in 1:length(pi_range)){
  
  ## OUTER LOOP OVER TRIALS
  reps<-10000
  trial_result<-vector(length=reps)
  trial_result_binary<-vector(length=reps)
  trial_result_ss<-vector(length=reps)
  p_hat<-vector(length=reps)
  inner_efficacy<-vector(length=reps)
  inner_futility<-vector(length=reps)
  inner_inconclusive<-vector(length=reps)
  inner_posterior_mean<-vector(length=reps)
  inner_coverage<-vector(length=reps)
  #freq_trial_result<-vector(length=reps)


  ## Step 3: Create inner loop for simulations of the specified trial design
  for (i in 1:reps){
    
    ## INNER LOOP OVER SAMPLE 
    efficacy<-0     # initialize
    futility<-0     # initialize
    inner_inconclusive[i]<-1 # initialize
    n<-1            # initializing sample size
    y1_before<-0    # initializing number of successes
    
    while (n<max_ss){
      
      # generate trial result
      y<-rbinom(1,size=1,prob=pi_range[j])
      # new number of successes
      y1<-sum(y1_before)+y
      # new number of failures
      y0<-n-y1
      # futility (even optimst would give up)
      futility<-pbeta(theta_futility,alpha_H+y1,beta_H+y0,lower.tail=TRUE)
      # efficacy (even pessimist would accept)
      efficacy<-pbeta(theta_L,alpha_L+y1,beta_L+y0,lower.tail=FALSE)
      ## for next loop
      y1_before<-y1

      if (n%%check_vector[k]==0 & futility>futility_threshold){
        inner_futility[i]<-1
        inner_inconclusive[i]<-0
        break
      }
      if (n%%check_vector[k]==0 & efficacy>efficacy_threshold){
        inner_efficacy[i]<-1
        inner_inconclusive[i]<-0
        break
      }
      
      n<-n+1
      }
    
    ## Posterior model probabilities, slide 133/1005 of 779 notes
    # posterior probability of data given models
    #post_L_D<-(beta(alpha_L,beta_L))^(-1)*beta(alpha_L+y1,beta_L+y0)
    #post_H_D<-(beta(alpha_H,beta_H))^(-1)*beta(alpha_H+y1,beta_H+y0)
    # posterior probability of models given data
    #post_L<-post_L_D*prior_L/(post_L_D*prior_L+post_H_D*prior_H)
    #post_H<-post_H_D*prior_H/(post_L_D*prior_L+post_H_D*prior_H)
    ## compute probability of event A given the models, slide 7/1005 of 779 notes
    
    
    post_L_eventA<-pbeta(theta_L,alpha_L+y1,beta_L+y0,lower.tail=FALSE)
    post_H_eventA<-pbeta(theta_L,alpha_H+y1,beta_H+y0,lower.tail=FALSE)
    # compute probability of event A
    #trial_result[i]<-sum(post_L*post_L_eventA,post_H*post_H_eventA)
    #trial_result_binary[i]=(trial_result[i]>=significant_threshold) # rejecting null hypothesis
    
    # inference posterior mixture proportions
    c1<-prior_L*beta(alpha_L+y1,beta_L+y0)/
      (prior_L*beta(alpha_L+y1,beta_L+y0)+prior_H*beta(alpha_H+y1,beta_H+y0))
    c2<-prior_H*beta(alpha_H+y1,beta_H+y0)/
      (prior_L*beta(alpha_L+y1,beta_L+y0)+prior_H*beta(alpha_H+y1,beta_H+y0))
    # posterior mean
    inner_posterior_mean[i]<-c1*(alpha_L+y1)/(alpha_L+beta_L+n)+c2*(alpha_H+y1)/(alpha_H+beta_H+n)
    
    x<-0
    for (l in 1:reps){
      p = runif(1)
      if(p<=c1){
        x[l] = rbeta(1,alpha_L+y1,beta_L+y0)
      }
      if(p>c1){
        x[l] = rbeta(1,alpha_H+y1,beta_H+y0)
      }
    }
    # coverage probability
    inner_coverage[i]<-(pi_range[j]>quantile(x,0.025) & pi_range[j]<quantile(x,0.975))
    
    trial_result_ss[i]<-n
    p_hat[i]<-y1/n
    #freq_trial_result[i]<-((sqrt(n)*(pi_range[j]-p_hat-delta_0)/sqrt(p_hat*(1-p_hat))<qnorm(alpha)) &
    #    (sqrt(n)*(pi_range[j]-p_hat+delta_0)/sqrt(p_hat*(1-p_hat))>qnorm(1-alpha)))
  }
  outer_trial_result_binary[k,j]<-mean(trial_result_binary) # sent to final results matrix
  outer_trial_result[k,j]<-mean(trial_result)
  #outer_freq_trial_result[j]<-mean(freq_trial_result)
  outer_trial_result_ss[k,j]<-mean(trial_result_ss)
  outer_p_hat[k,j]<-mean(p_hat)
  outer_futility[k,j]<-mean(inner_futility)
  outer_efficacy[k,j]<-mean(inner_efficacy)
  outer_inconclusive[k,j]<-mean(inner_inconclusive)
  outer_posterior_mean[k,j]<-mean(inner_posterior_mean)
  outer_coverage[k,j]<-mean(inner_coverage)
}
}

par(ask=FALSE)
plot(pi_range,outer_futility,type='l',ylim=c(0,1))
lines(pi_range,outer_efficacy)
lines(pi_range,outer_inconclusive)
outer_trial_result_ss
outer_p_hat
outer_posterior_mean
outer_coverage


.outer_trial_result_binary # probability of rejecting null hypothesis (alpha at theta_L, 1-beta at theta_H)
outer_trial_result
outer_trial_result_ss

## Plot of power curves
par(ask=TRUE)
plot(pi_range,outer_trial_result_binary[1,],type='l',
     xlab="Response Probability",ylab="Probability of Rejection",main="Power Curves")
lines(pi_range,outer_trial_result_binary[2,],type='l')
lines(pi_range,outer_trial_result_binary[3,],type='l')
lines(pi_range,outer_trial_result_binary[4,],type='l')
legend(x=theta_H*(2/3),y=0.8,
       title="# enrolled",
       legend=check_vector,
       col=check_vector,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=TRUE)
plot(pi_range,outer_trial_result_ss[4,],type='l',
     xlab="Response Probability",ylab="Expected Sample Size",main="Sample Size",
     ylim=c(max_ss/2,max_ss*1.1))
axis(2,at=max_ss,labels="Max")
points(pi_range,outer_trial_result_ss[3,],type='l',lty=2)
points(pi_range,outer_trial_result_ss[2,],type='l',lty=3)
points(pi_range,outer_trial_result_ss[1,],type='l',lty=4)
legend(x=theta_L*(3/2),max_ss*(4/5),
       title="# enrolled",
       legend=check_vector,
       col=check_vector,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=TRUE)
plot(pi_range,abs(outer_p_hat[4,]-pi_range),type='l',
     xlab="Response Probability",ylab="Bias",main="Bias",
     ylim=c(0,0.1))
points(pi_range,abs(outer_p_hat[3,]-pi_range),type='l',lty=2)
points(pi_range,abs(outer_p_hat[2,]-pi_range),type='l',lty=3)
points(pi_range,abs(outer_p_hat[1,]-pi_range),type='l',lty=4)
legend(x=theta_H*(6/7),0.1,
       title="# enrolled",
       legend=check_vector,
       col=check_vector,lty=c(4,3,2,1),
       cex=0.8,text.font=3)
par(ask=FALSE)

#}

#example1()
