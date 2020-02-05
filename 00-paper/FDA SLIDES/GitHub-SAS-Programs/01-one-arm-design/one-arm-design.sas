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
%include "&macPath.\violin_se.sas";
%include "&macPath.\power.sas";

%let outPath = P:\Projects\FDA-IPA\Bayesian-Sequential-Monitoring\Presentations\2018-12-14\figures;

libname out "&prgPath\data";

%let PiNull = 0.20;
%let PiAlt  = 0.40;

%let mu0    = 0.20;
%let phi0   = 14.0;

%let mu1    = 0.40;
%let phi1   = 14.0;

%let nby    =    2;


ods select none;
ods output Output = PowerAnalysis;
proc power;
   onesamplefreq test=adjz method=normal
      nullproportion = &piNull.
      proportion = &piAlt.
      sides = u
      ntotal = .
      power = 0.90;
run;
ods select all;

proc SQL noprint;
 select 76 into :MaxN from PowerAnalysis;
quit;
%put &=MaxN;

%let FU_DIST    = rand('Normal',4,0.25);
%let ER_DIST    = rand('gamma',1,0.50); 
%let DSOUT      = VIOLIN;
%sim_se(seed=3,nSims=1,plots=NONE);

%violin_se(pi=0.20,label=NULL);
%violin_se(pi=0.40,label=ALT);


%let FU_DIST    = rand('Normal',4,0.25);
%let ER_DIST    = rand('gamma',1,0.50); 
%let DSOUT      = slow;
%sim_se(seed=7861,nSims=50000);
%plotA;
%plotB;

proc SQL noprint;
 select 50 into :MaxN from PowerAnalysis;
quit;
%put &=MaxN;


%let FU_DIST    = rand('Normal',4,0.25);
%let ER_DIST    = rand('gamma',1,0.125); 
%let DSOUT      = fast;
%sim_se(seed=2165,nSims=50000);
%plotA;
%plotB;

/**/
/**/
/*%let FU_DIST    = rand('Normal',4,0.25);*/
/*%let ER_DIST    = rand('gamma',1,0.10); */
/*%let DSOUT      = veryFast;*/
/*%sim_se(seed=4264,nSims=50000);*/
/*%plotA;*/
/*%plotB;*/












