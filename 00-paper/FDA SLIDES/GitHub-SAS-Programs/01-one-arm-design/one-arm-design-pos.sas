/*****************************************************************************
* Project           : FDA Consulting
*
* Program name      : one-arm-design-pos.sas
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

%let root = P:\Projects\FDA-IPA\Bayesian-Sequential-Monitoring;
%let prgPath = &root.\programs\01-one-arm-design;
%let macPath = &prgPath\macros;

%include "&macPath.\sim_pos.sas";
%include "&macPath.\violin_pos.sas";
%include "&macPath.\power_pos.sas";

%let outPath = P:\Projects\FDA-IPA\Bayesian-Sequential-Monitoring\Presentations\2018-12-14\figures;
libname out "&prgPath\data";

ods noresults;
ods listing close;

%let PiNull = 0.20;
%let PiAlt  = 0.40;


%let mu0  = 0.20;
%let phi0 = 15.0;

%let mu1  = 0.40;
%let phi1 = 14.0;


%let nby    =    2;

proc SQL noprint;
 select 76 into :MaxN from sashelp.class;
quit;
%put &=MaxN;

%let FU_DIST    = rand('Normal',4,0.25);
%let ER_DIST    = rand('gamma',1,0.50); 
%let DSOUT      = slow;
%sim_pos(seed=2,nSims=1,plots=NONE);
%violin_pos(pi=0.20,label=pos-NULL);




%let FU_DIST    = rand('Normal',4,0.25);
%let ER_DIST    = rand('gamma',1,0.50); 
%let DSOUT      = slow;
%sim_pos(seed=7861,nSims=50000);
%plotA;
%plotB;






/**/
/*data check;*/
/**/
/*   yObs = 22;*/
/*   nObs = 50;*/
/**/
/*   %macro pos(pref=int,nVar=nObs,yVar=yObs);*/
/*   nFill = &maxN. - &nVar.;*/
/**/
/*    &pref._postProbOpt=0;*/
/*	do x = 0 to nFill;*/
/*	   logPMF    = lcomb(nFill, x) +*/
/*	               logbeta(x + &yVar. + &mu0.*&phi0., nFill - x + &nVar.-&yVar. + (1-&mu0.)*&phi0.) -*/
/*	               logbeta(&yVar. + &mu0.*&phi0., &nVar.-&yVar. + (1-&mu0.)*&phi0.);*/
/*	   PMF       = exp(logPMF);        */
/*	   INTEGRAND = (sdf('beta',&PiNull.,x+&yVar. + &mu0.*&phi0. ,nFill+&nVar.-&yVar.-x + (1-&mu0.)*&phi0.)>=0.95);*/
/*	   &pref._postProbOpt + INTEGRAND * PMF; */
/*    output; */
/*	end;*/
/*    &pref._fut           = (&pref._postProbOpt   <= 0.05)*(nObs>=20);*/
/*   %mend pos;*/
/*   %pos(pref=int,nVar=nObs,yVar=yObs);*/
/*run;*/
