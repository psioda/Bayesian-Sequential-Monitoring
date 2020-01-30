#########################################
#### Figure 5, Risk Diff Prior Plots ####
#########################################

library(rmutil)
library(lattice)
library(pracma)
source("../args_model.R")
source("../code_functions.R")

fcn_prior_placebo()
skpt_tail_area()
enth_tail_area()

## SKEPTICAL PRIORS ##
prior.skpt<-function(x,y){
  exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
    exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)}
nc.skpt<-tryCatch(integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
             error=function(e) integral2(prior.skpt,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
prior.skpt.nc<-function(x,y){
  exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
    exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt}

## ENTHUASTIC PRIORS ##
prior.enth<-function(x,y){
  exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
    exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)}
nc.enth<-tryCatch(integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,singular=T)[[1]],
             error=function(e) integral2(prior.enth,xmin=0,xmax=1,ymin=0,ymax=1,abstol=1E-6)[[1]])
prior.enth.nc<-function(x,y){
  exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
    exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/nc.enth}

x.len<-101
grid1d <- seq(0, 1, length= x.len)
skpt.x<-0
skpt.y<-0
enth.x<-0
enth.y<-0

for (i in 1:length(grid1d)){
  ## SKEPTICAL PRIORS ##
  x<-grid1d[i]
  marginal.x<-function(y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt}
  skpt.x[i]<-integrate(marginal.x,lower=0,upper=1)[[1]]
  y<-grid1d[i]
  marginal.y<-function(x){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.skpt)/sigma0.skpt)^lambda0.skpt)/nc.skpt 
  }
  skpt.y[i]<-integrate(marginal.y,lower=0,upper=1)[[1]]
  ## ENTHUASTIC PRIORS ##
  x<-grid1d[i]
  marginal.x<-function(y){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/nc.enth}
  enth.x[i]<-integrate(marginal.x,lower=0,upper=1)[[1]]
  y<-grid1d[i]
  marginal.y<-function(x){
    exp(-(abs(x-mu)/sigma0.placebo)^lambda0.placebo)*
      exp(-(abs((y-x)-delta.enth)/sigma0.enth)^lambda0.enth)/nc.enth 
  }
  enth.y[i]<-integrate(marginal.y,lower=0,upper=1)[[1]]
}