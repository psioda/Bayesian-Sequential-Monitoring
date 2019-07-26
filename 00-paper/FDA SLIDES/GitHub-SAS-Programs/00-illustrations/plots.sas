/*****************************************************************************
* Project           : FDA Consulting
*
* Program name      : one-arm-design
*
* Author            : Matthew A. Psioda
*
* Date created      : 2018-12-07
*
* Purpose           : This program is the solution to final exam part 1
*
* Revision History  :
*
* Date          Author   Ref (#)  Revision
* 
*
*
* searchable reference phrase: *** [#] ***;
******************************************************************************/
*ds noresults;
ods listing close;

%let root = P:\Projects\FDA-IPA\Bayesian-Sequential-Monitoring;
%let prgPath = &root.\programs\00-illustrations;


%let outPath = P:\Projects\FDA-IPA\Bayesian-Sequential-Monitoring\Presentations\2018-12-14\figures;




%let PiNull = 0.20;
%let PiAlt  = 0.40;

%let mu0  = 0.20;
%let phi0 = 14.0;

%let mu1  = 0.40;
%let phi1 = 14.0;

ods html file = "&outPath.\freqPower.html";
ods trace on;
ods output Output = PowerAnalysis;
proc power;
   onesamplefreq test=adjz method=normal
      nullproportion = &piNull.
      proportion = &piAlt.
      sides = u
      ntotal = .
      power = 0.90;
run;
ods trace off;
ods html close;

data plot;
 do pi = 0.001 to 0.999 by 0.001;
   if pi <= &PiAlt. then pL=0; else pL = .;
   if pi <= &PiAlt. then pH=.; else pH = 0;
   Prior = 'Skeptical   '; pdf = pdf('beta',pi,&mu0.*&phi0.,(1-&mu0.)*&phi0.); output;

   if pi <= &PiNull. then pL=0; else pL = .;
   if pi <= &PiNull. then pH=.; else pH = 0;

   Prior = 'Enthusiastic'; pdf = pdf('beta',pi,&mu1.*&phi1.,(1-&mu1.)*&phi1.); output;
 end;
run;

data tails;
 pSkept = sdf('beta',&PiAlt,&mu0.*&phi0.,(1-&mu0.)*&phi0.);   call symput('pSkept',strip(put(pSkept,7.3)));
 pOpt   = cdf('beta',&PiNull, &mu1.*&phi1.,(1-&mu1.)*&phi1.); call symput('pOpt',strip(put(pOpt,7.3)));
run;




ods noproctitle;
title;
options topmargin=0.1in bottommargin=0.1in rightmargin=0.1in leftmargin=0.1in nodate nonumber;
option papersize=("5.4in","6.4in") nodate nonumber;

ods pdf file = "&outPath.\skeptical.pdf" startpage=no;
 
ods graphics / height=5in width=5in noborder;

proc odstext; p "Skeptical Prior" / style=[font=(Arial) fontsize=9pt just=c]; run;
proc sgplot data = plot noautolegend;
 where Prior = 'Skeptical';
 band x=pi lower=pL upper=pdf / fillattrs=(color=blue);
 band x=pi lower=pH upper=pdf / fillattrs=(color=Red);
 yaxis label = 'Density Value';
 xaxis label = 'Response Probability' display=( noticks novalues nolabel);
run;
quit;


ods pdf close;


ods noproctitle;
title;
options topmargin=0.1in bottommargin=0.1in rightmargin=0.1in leftmargin=0.1in nodate nonumber;
option papersize=("5.4in","6.4in") nodate nonumber;

ods pdf file = "&outPath.\enthusiastic.pdf" startpage=no;
 
ods graphics / height=5in width=5in noborder;

proc odstext; p "Enthusiastic Prior" / style=[font=(Arial) fontsize=9pt just=c]; run;
proc sgplot data = plot noautolegend;
 where Prior = 'Enthusiastic';
 band x=pi lower=pL upper=pdf / fillattrs=(color=red);
 band x=pi lower=pH upper=pdf / fillattrs=(color=blue);
 yaxis label = 'Density Value';
 xaxis label = 'Response Probability' display=( noticks novalues nolabel);
run;
quit;


ods pdf close;
