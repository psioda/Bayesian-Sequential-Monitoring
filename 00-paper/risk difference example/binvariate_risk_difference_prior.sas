%let sigma0 = 3.5;
%let mu0    = 0.50;
%let alpha0 = 20;

%let sigma1 = 7;
%let alpha1 = 15;
%let delta1 = 0.00;

proc IML ;

 start den01(z);
  L0 = -(  (&sigma0.*( z[,2]-&mu0.            ))##2     )##&alpha0.;
  L1 = -(  (&sigma1.*( z[,1]-(z[,2]+&delta1.) ))##2     )##&alpha1.;
  return(  exp( L0 + L1 ) );
 finish;

 start den0(x);
  L0 = -(   (&sigma0.*(x-&mu0.))##2 )##&alpha0.;
  return(  exp( L0) );
 finish;

 interval = 0||1;

 * normalizing constant for den0;
 nc = .;
 call quad(nc,'den0',interval);
 print nc;

 * normalizing constant for den01;
 by = 0.025; * default 0.00025 for valid integration;
 x = t(do(0.0000,1.0000,by));
 y = x;
 z = expandGrid(y,x);
 f = den01(z);
 nc01 = sum(f*by*by);
 print nc01;
 call symput('nc01',strip(char(nc01)));

 create grid from z[colname={"x" "y"}]; 
 append from z; 
 close grid;

 create height var{f}; 
 append; 
 close height;

 call series(y,den0(y)/nc);

quit;

data plot;
merge grid height;
run;


goptions reset=all border cback=white htitle=12pt; 
proc g3d data=plot; plot x*y=f; run; quit;

data x;
run;

proc mcmc data = x nmc=50000 nbi=1000 mintune=10 ntu=1000 postout = post;

 array pi[2] pi0 pi1 (0.30 0.30); * robust to choice of initial value;
 parms pi1 pi0 / slice;

 beginprior;
  L0 = -(  (&sigma0.*( pi0-&mu0.            ))**2     )**&alpha0.;
  L1 = -(  (&sigma1.*( pi1-(pi0+&delta1.)   ))**2     )**&alpha1.;
  logPrior = L0 + L1 - &nc01.;
 endprior;

 prior pi1 ~ general(0,lower=0,upper=1);  *same as uniform(0,1);
 prior pi0 ~ general(0,lower=0,upper=1);  *same as uniform(0,1);

 model general(logPrior);

run;
quit;

data post;
 set post;
 byValue = 0.025;
 do j = 0 to 1 by byValue;
  if j <= pi0 < j + byValue then do;
   category = j + byValue*0.5;
   j = 1000;
  end;
 end;
run;

proc sgplot data = post;
 histogram pi0 / nbins=100 binstart=0;
 xaxis values = (0 to 1 by 0.10);
run;

proc sgplot data = post;
 histogram pi1 / nbins=100 binstart=0;
 xaxis values = (0 to 1 by 0.10);
run;

proc sgplot data = post;
 vbox pi1 / category=category;
run;
