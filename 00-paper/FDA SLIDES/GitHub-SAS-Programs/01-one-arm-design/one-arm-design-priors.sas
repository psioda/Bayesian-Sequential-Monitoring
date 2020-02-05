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
ods noresults;
ods listing close;

%let root = P:\Projects\FDA-IPA\Bayesian-Sequential-Monitoring;
%let prgPath = &root.\programs\01-one-arm-design;
%let macPath = &prgPath\macros;

%include "&macPath.\sim_se.sas";

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
option papersize=("8.8in","5.4in") nodate nonumber;

ods pdf file = "&outPath.\one-arm-priors.pdf";
 
ods layout start height=5.20in width=8.60in; 

ods region x=0.1in y=0.1in width=4.0in height=5.0in; 
ods graphics / height=4in width=4in;
proc odstext; p "Skeptical Beta Prior" / style=[font=(Arial) fontsize=9pt just=c]; run;
proc sgplot data = plot noautolegend;
 where Prior = 'Skeptical';
 band x=pi lower=pL upper=pdf / fillattrs=(color=blue);
 band x=pi lower=pH upper=pdf / fillattrs=(color=Red);
 yaxis label = 'Density Value';
 xaxis label = 'Response Probability';
 inset "P((*ESC*){unicode theta}>&PiAlt.)=&pSkept." "Approximately N=&phi0 Subjects" "E[(*ESC*){unicode theta}]=&PiNull." / position=topright;
run;
quit;

ods region x=4.2in y=0.1in width=4.0in height=5.0in; 
ods graphics / height=4in width=4in;
proc odstext; p "Enthusiastic Beta Prior" / style=[font=(Arial) fontsize=9pt just=c]; run;
proc sgplot  data = plot  noautolegend;
 where Prior = 'Enthusiastic';
 band x=pi lower=pL upper=pdf / fillattrs=(color=Red);
 band x=pi lower=pH upper=pdf / fillattrs=(color=Blue);
 yaxis display=(nolabel);
 xaxis label = 'Response Probability';
 inset "P((*ESC*){unicode theta}<&PiNull.)=&pOpt." "Approximately N=&phi1 Subjects" "E[(*ESC*){unicode theta}]=&PiAlt." / position=topright;
run;
quit;

ods layout end;

ods pdf close;


data plot;
 do pi = 0.001 to 0.999 by 0.001;
   p = 0;
   Prior = 'Mixture'; pdf = 0.5*pdf('beta',pi,&mu0.*&phi0.,(1-&mu0.)*&phi0.)
                          + 0.5*pdf('beta',pi,&mu1.*&phi1.,(1-&mu1.)*&phi1.); 

   if pi >= 0.20 and pi <= 0.40 then pdf2 = pdf;
   else pdf2 = .;
   output;
 end;
run;

ods noproctitle;
title;
options topmargin=0.1in bottommargin=0.1in rightmargin=0.1in leftmargin=0.1in nodate nonumber;
option papersize=("5in","5in") nodate nonumber;

ods pdf file = "&outPath.\one-arm-priors-mixture.pdf" startpage=no;
proc odstext; p "50/50 Mixture of Skeptical and Enthusiastic Priors" / style=[font=(Arial) fontsize=9pt just=c]; run;
ods graphics / height=4in width=4in;
proc sgplot data = plot noautolegend;
 band x=pi lower=p upper=pdf  / fillattrs=(color=blue transparency=0.2);
 band x=pi lower=p upper=pdf2 / fillattrs=(color=blue transparency=0.6)  y2axis;

 yaxis label = 'Density Value';
 y2axis display=(nolabel);
 xaxis label = 'Response Probability';
run;
quit;

ods pdf close;
