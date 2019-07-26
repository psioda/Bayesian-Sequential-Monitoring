%macro violin_pos(pi=0.40,label=ALT);
data example_plot;
 set example;
 by pi;
 where round(pi,1e-5)=&pi.;
 if final = . then final = 0;

 
 
 
   lower_cr = 0.99999;
   tail  = 1;
   do while(tail>0.025);
	 tail = pw0*cdf('beta',lower_cr,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.)
		  + pw1*cdf('beta',lower_cr,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);
		  
		  if tail >= 0.025 then lower_cr = lower_cr - 1e-5;
   end;

   upper_cr = 0.00001;
   tail  = 1;
   do while(tail>0.025);
	 tail = pw0*sdf('beta',upper_cr,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.)
		  + pw1*sdf('beta',upper_cr,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);
		  
		  if tail >= 0.025 then upper_cr = upper_cr + 1e-4;
   end;		
 
 
 
 do p = 0.001 to 0.800 by 0.001;
  Prior = 'Skeptical  '; density = pdf('beta',p,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.); output;
  Prior = 'Enthusiatic'; density = pdf('beta',p,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.); output;
 end;

 if last.pi then do;
  n = 0;
  nObs = 0;
  yObs = 0;
  nMis = 0;
  final = 0;

  int_mean = 0.5*&mu0. + 0.5*&mu1.;
  
  
   lower_cr = 0.99999;
   tail  = 1;
   do while(tail>0.025);
	 tail = pw0*cdf('beta',lower_cr,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.)
		  + pw1*cdf('beta',lower_cr,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);
		  
		  if tail >= 0.025 then lower_cr = lower_cr - 1e-5;
   end;

   upper_cr = 0.00001;
   tail  = 1;
   do while(tail>0.025);
	 tail = pw0*sdf('beta',upper_cr,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.)
		  + pw1*sdf('beta',upper_cr,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);
		  
		  if tail >= 0.025 then upper_cr = upper_cr + 1e-4;
   end;		
  
  
  
 
	 do p = 0.001 to 0.800 by 0.001;
	  Prior = 'Skeptical  '; density = pdf('beta',p,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.); output;
	  Prior = 'Enthusiatic'; density = pdf('beta',p,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.); output;
	 end;
 end;
run; proc sort; by prior final n nObs yObs nMis int_mean; run;



data example_plot;
 set example_plot;
 by prior final n nObs yObs nMis int_mean;
  if not first.nMis then do;
    lower_cr = .;
	upper_cr = .;  
  end;
run;

proc SQL noprint;
 select max(density) into : maxDen from example_plot;
 select count( distinct nObs) into : maxTPT from example_plot;

quit;

%put &maxTPT.;

data example_plot_dens2;
  merge example_plot(where=(Prior='Skeptical'  ) rename=(p=p1 density=density1))
        example_plot(where=(Prior='Enthusiatic') rename=(p=p2 density=density2));
  by n nobs;

  if first.nobs then index+1;
  if not first.nobs then do;
   int_mean = .;
   n = .;
   nobs = .;
   nmis = .;
   yobs = .;
   int_postProbOpt = .;
  end;

        density2 =  density2 / (2.1*&maxDen);
		density1 =  density1 / (2.1*&maxDen);


	    zero = index;

		array densityA[&maxTPT.];
		array densityB[&maxTPT.];

		densityA[index] = index - density1;
		densityB[index] = index + density2;


	if final = 0 then output;
	if final = 1 then do;
     p1b = p1; p1 = .;
	 p2b = p2; p2 = .;
	 output;
	end;
run; 

ods noproctitle;
title;
options topmargin=0.1in bottommargin=0.1in rightmargin=0.1in leftmargin=0.1in nodate nonumber;
option papersize=("8.2in","4.2in") nodate nonumber;
ods graphics / height=4in width=8in noborder;
ods pdf file = "&outPath.\one-arm-design-example-path-&label..pdf";


proc sgplot data=example_plot_dens2 nocycleattrs;
/*  sgplot  nObs yObs nMis / rows=1 onepanel */
/*    novarname noborder colheaderpos=bottom;*/

 refline 1.5 10.5 %sysevalf(&maxTPT.-0.5) / axis=x lineattrs=(pattern=2 color=gray);
 %macro A;
  %do j = 1 %to &maxTPT.;
  band y=p1 upper=zero lower=densityA&j.    / fill fillattrs=(color=red  transparency=0.5)  /*outline lineattrs=(color=black pattern=1)*/ name="A" ;
  band y=p2 upper=densityB&j. lower=zero    / fill fillattrs=(color=blue transparency=0.5) /*outline lineattrs=(color=black pattern=1)*/ name="B";

  band y=p1b upper=zero lower=densityA&j.   / fill fillattrs=(color=Darkred transparency=0.5)  /*outline lineattrs=(color=black pattern=1)*/ name = "C";
  band y=p2b upper=densityB&j. lower=zero   /  fill fillattrs=(color=Darkblue transparency=0.5) /*outline lineattrs=(color=black pattern=1)*/ name = "D";
  %end;
 %mend A;
 %A;
  series x=index y=int_mean                  / markers markerattrs=(symbol=diamondFilled color=purple) lineattrs=(color=purple thickness=1.5); 
  highlow x=index low=lower_cr high=upper_cr / lineattrs=(color=purple thickness=1.5) highcap=serif lowcap=serif name = "E" legendLabel='Mixture Prior 95% Credible Interval';

  format int_postProbOpt 6.3;
  xaxistable int_postProbOpt / label='PoS';
  
/*  xaxistable n    / label="n (ovr.)";*/
  xaxistable nobs / label='n (obs.)';
  xaxistable yobs / label='y (obs.)';
  xaxistable nmis / label='y (mis.)';

  yaxis  grid label = 'Probability of Response';
  label p1 = 'Skeptical Prior Posterior'
        p2 = 'Enthusiastic Prior Posterior';
		keylegend "A" "B" "E";
  xaxis values=(0.5 to %sysevalf(&maxTPT.+0.5) by 0.5) display=none;
run;

ods pdf close;
%mend violin_pos;
