libname outdata "D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/SAS";


ods noproctitle;

ods html close;
ods html newfile=none;
ods html;

proc IML;

mu0 = 0.4;
sigma0 = 0.1;
lambda0 = 2;

n = 5;
y = 2;

scaleValue = 1e-6;

start den(x) global(skip_data,y,n,mu0,sigma0,lambda0,log_nc);
  if skip_data = 0 then      return exp( y*log(x) + (n-y)*log(1-x) - (  (abs(x-mu0)/sigma0)##lambda0 ) - log_nc);
  else if skip_data = 1 then return exp(  - 0.5*(  (abs(x-mu0)/sigma0)##lambda0 ) - log_nc);
finish;

/** compute prior normalizing constant **/
skip_data = 1;
log_nc    = 0;
interval  = 1e-6||0.999999;
call quad(computed_log_nc,"den",interval) eps=1e-11 peak=mu0;
log_nc = log(computed_log_nc);

/** plot prior **/
x = t(do(0.0001,0.9999,0.001));
p = den(x);
title "Prior";
run series(x,p);

skip_data = 0;
interval  = 1e-6||0.999999;
call quad(marg_like,"den",interval) eps=1e-11;
log_nc = log_nc + log(marg_like);

p = den(x);
title "Posterior";
call series(x,p);

interval=0.20||0.999999;
call quad(post_prob,"den",interval) eps=1e-11;

title;
print post_prob[l="Posterior Probability"];

quit;

data person;
   infile datalines delimiter=','; 
   input y;
   datalines;                      
0
1
0
1
0
;
run;
proc print data=person; run;

proc mcmc data=person seed=1 nbi=5000 nmc=100000 outpost=regOut;
parms p;
prior p ~ beta(1,1);
model y~binary(p);
run;

# 0.4283
# by hand, posterior is beta(3,4), with expected value 3/7=0.42857142857
