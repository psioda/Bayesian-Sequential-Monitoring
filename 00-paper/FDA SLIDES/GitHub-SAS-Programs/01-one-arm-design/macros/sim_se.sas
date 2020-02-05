%macro sim_se(seed=1,nSims=10000,plots=ALL);

data sim(keep=true_pi sim n yObs nObs nMis yFin nFin int_postProbSkept int_eff int_postProbOpt int_fut int_mean int_lower_cr int_upper_cr int_cover
                                                     fin_postProbSkept fin_eff fin_postProbOpt fin_fut fin_mean fin_lower_cr fin_upper_cr fin_cover)
 example(keep=pi final n yObs nObs nMis int_postProbSkept int_postProbOpt int_mean int_eff int_fut pw1 pw0);
 call streaminit(&seed.);
 Ntot = &MaxN.;

 array r[&maxN.];
 array e[&maxN.];
 array y[&MaxN.];
 array ord[&maxN.];
 


	 do pi = 0.15 to 0.45 by 0.05;
     do sim = 1 to &nSims.;

  
	     ** fill out all the hypothetical data;
	     r[1] = 0;
		 do j = 1 to dim(r);
		  if j>1 then r[j] = r[j-1] + &ER_DIST.;
		  e[j]   = r[j] + &FU_DIST.;
		  ord[j] = e[j];
		  y[j]   = rand('bernoulli',pi);
		 end;
		 call sortn(of ord[*]);


 
		 ** perform the sequential analysis;
		 stop = 0;
         n = 0;
		 do while(stop=0);
          n = n + &nBy.;
		  yObs  = 0;
		  nObs  = 0;
		  nMis  = 0;
		  yFin  = 0;
		  nFin  = 0;
		  final = 0;

		  eTime = ord[n];
		  do j = 1 to dim(r);
		        if e[j] <= eTime then do; yObs + y[j]; nObs + 1; end;
		   else if r[j] <= eTime then do; nMis + 1; end;
		  end;

		   int_postProbSkept = sdf('beta',&PiNull.,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.);
		   int_postProbOpt   = cdf('beta',0.5*&PiAlt. + 0.5*&PiNull., yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);

		   int_eff = (int_postProbSkept >= 0.95);
		   int_fut = (int_postProbOpt   >= 0.85);


		   int_mean0 = (yObs + &mu0.*&phi0.) / (yObs + &mu0.*&phi0.  + nObs-yObs + (1-&mu0.)*&phi0.);
		   int_mean1 = (yObs + &mu1.*&phi1.) / (yObs + &mu1.*&phi1.  + nObs-yObs + (1-&mu1.)*&phi1.);

		   w0 = lgamma(       &phi0.)  - lgamma(&mu0.*&phi0.)        - lgamma(            (1-&mu0.)*&phi0.)
		      - lgamma(nObs + &phi0.)  + lgamma(yObs + &mu0.*&phi0.) + lgamma(nObs-yObs + (1-&mu0.)*&phi0.);

		   w1 = lgamma(       &phi1.)  - lgamma(&mu1.*&phi1.)        - lgamma(            (1-&mu1.)*&phi1.)
		      - lgamma(nObs + &phi1.)  + lgamma(yObs + &mu1.*&phi1.) + lgamma(nObs-yObs + (1-&mu1.)*&phi1.);

		   ww0 = exp(w0 - max(w0,w1));
           ww1 = exp(w1 - max(w0,w1));

		   pw0 = ww0 / (ww0+ww1);
		   pw1 = ww1 / (ww0+ww1);

           int_mean = pw0*int_mean0 + pw1*int_mean1;




		   if int_eff    then stop = 1;
		   if int_fut    then stop = 1;
		   if n = &MaxN. then stop = 1;

           if sim = 1 and round(pi,1e-5) in (0.20,0.40) then do;
		      output example;
           end;
		   
		   
		   if stop = 1 then do;
			   lower_cr = 0.60;
			   tail  = 1;
			   do while(tail>0.025);
				 tail = pw0*cdf('beta',lower_cr,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.)
					  + pw1*cdf('beta',lower_cr,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);
					  
					  if tail >= 0.025 then lower_cr = lower_cr - 1e-3;
			   end;

			   upper_cr = 0.10;
			   tail  = 1;
			   do while(tail>0.025);
				 tail = pw0*sdf('beta',upper_cr,yObs + &mu0.*&phi0. ,nObs-yObs + (1-&mu0.)*&phi0.)
					  + pw1*sdf('beta',upper_cr,yObs + &mu1.*&phi1. ,nObs-yObs + (1-&mu1.)*&phi1.);
					  
					  if tail >= 0.025 then upper_cr = upper_cr + 1e-3;
			   end;	
				int_lower_cr = lower_cr;
				int_upper_cr = upper_cr;
				
				if lower_cr <= pi <= upper_cr then int_cover=1;
				else int_cover = 0;
		   end;
		   
     	 end;

         ** perform analysis when all outcomes observed;
		  yFin  = 0;
		  nFin  = 0;
		  do j = 1 to dim(r);
            if r[j] <= eTime then do; yFin + y[j]; nFin + 1; end;
		  end;

		   fin_postProbSkept = sdf('beta',&PiNull.,yFin + &mu0.*&phi0. ,nFin-yFin + (1-&mu0.)*&phi0.);
		   fin_postProbOpt   = cdf('beta',0.5*&PiAlt. + 0.5*&PiNull., yFin + &mu1.*&phi1. ,nFin-yFin + (1-&mu1.)*&phi1.);

		   fin_eff = (fin_postProbSkept >= 0.95);
		   fin_fut = (fin_postProbOpt   >= 0.85);

		   fin_mean0 = (yFin + &mu0.*&phi0.) / (yFin + &mu0.*&phi0.  + nFin-yFin + (1-&mu0.)*&phi0.);
		   fin_mean1 = (yFin + &mu1.*&phi1.) / (yFin + &mu1.*&phi1.  + nFin-yFin + (1-&mu1.)*&phi1.);

		   w0 = lgamma(       &phi0.)  - lgamma(&mu0.*&phi0.)        - lgamma(            (1-&mu0.)*&phi0.)
		      - lgamma(nFin + &phi0.)  + lgamma(yFin + &mu0.*&phi0.) + lgamma(nFin-yFin + (1-&mu0.)*&phi0.);

		   w1 = lgamma(       &phi1.)  - lgamma(&mu1.*&phi1.)        - lgamma(            (1-&mu1.)*&phi1.)
		      - lgamma(nFin + &phi1.)  + lgamma(yFin + &mu1.*&phi1.) + lgamma(nFin-yFin + (1-&mu1.)*&phi1.);

		   ww0 = exp(w0 - max(w0,w1));
           ww1 = exp(w1 - max(w0,w1));

		   pw0 = ww0 / (ww0+ww1);
		   pw1 = ww1 / (ww0+ww1);

           fin_mean = pw0*fin_mean0 + pw1*fin_mean1;

		   
		   

			   lower_cr = 0.60;
			   tail  = 1;
			   do while(tail>0.025);
				 tail = pw0*cdf('beta',lower_cr,yFin + &mu0.*&phi0. ,nFin-yFin + (1-&mu0.)*&phi0.)
					  + pw1*cdf('beta',lower_cr,yFin + &mu1.*&phi1. ,nFin-yFin + (1-&mu1.)*&phi1.);
					  
					  if tail >= 0.025 then lower_cr = lower_cr - 1e-3;
			   end;

			   upper_cr = 0.10;
			   tail  = 1;
			   do while(tail>0.025);
				 tail = pw0*sdf('beta',upper_cr,yFin + &mu0.*&phi0. ,nFin-yFin + (1-&mu0.)*&phi0.)
					  + pw1*sdf('beta',upper_cr,yFin + &mu1.*&phi1. ,nFin-yFin + (1-&mu1.)*&phi1.);
					  
					  if tail >= 0.025 then upper_cr = upper_cr + 1e-3;
			   end;	
			   fin_lower_cr = lower_cr;
			   fin_upper_cr = upper_cr;
			   
  				if lower_cr <= pi <= upper_cr then fin_cover=1;
				else fin_cover = 0; 
		   
		   
		   
		   
		   
		   
		   
           if sim = 1 and round(pi,1e-5) in (0.20,0.40) then do;
             final = 1;
			 yObs = yFin;
			 nObs = nFin;
			 nMis = 0;
             int_postProbSkept = fin_postProbSkept;
			 int_postProbOpt   = fin_postProbOpt;
			 int_mean          = fin_mean;
             output example;
		   end;

		 true_pi = pi;
         output sim;
	  end;
      end;


run;

data out.example_&dsout.;
 set example;
run;


data Sim out.sim_se_&dsout.;
 set sim;

 agree_eff = (int_eff=fin_eff);
 agree_fut = (int_fut=fin_fut);

      length lev $20;

      if fin_postProbSkept >  0.95 then Lev = "(0.95,1.00)";
 else if fin_postProbSkept >  0.90 then Lev = "(0.90,0.95]";
 else if fin_postProbSkept >  0.80 then Lev = "(0.80,0.90]";
 else if fin_postProbSkept >  0.70 then Lev = "(0.70,0.80]";
 else if fin_postProbSkept >  0.50 then Lev = "(0.50,0.70]";
 else                                   Lev = "(0.00,0.50]";
run;

proc means data = sim noprint nway;
 class true_pi;
 var n yObs nObs nMis yFin nFin int_postProbSkept int_eff fin_postProbSkept fin_eff  
                                int_postProbOpt   int_fut fin_postProbOpt   fin_fut
                                agree_eff agree_fut int_mean fin_mean int_cover fin_cover;
 output out = summary(drop=_:)
  mean=n yObs nObs nMis yFin nFin int_postProbSkept int_eff fin_postProbSkept fin_eff  
                                int_postProbOpt   int_fut fin_postProbOpt   fin_fut
                                agree_eff agree_fut int_mean fin_mean int_cover fin_cover;
run;

proc means data = sim noprint nway;
 class true_pi;
 where agree_eff = 0 and int_eff=1;
 var fin_postProbSkept ;
 output out = summary_disagree1(drop=_:)
  mean = diagree_pp p10 = diagree_pp_L p90= diagree_pp_U;
run;

/*proc means data = sim noprint nway;*/
/* class true_pi;*/
/* where agree_fut = 0 and fin_fut=1;*/
/* var fin_postProbOpt ;*/
/* output out = summary_disagree2(drop=_:)*/
/*  mean = diagree_pp2 p10 = diagree_pp2_L p90= diagree_pp2_U;*/
/*run;*/



proc freq data = sim noprint;
 by true_pi ;
 table agree_eff / out = summary_eff_agree(drop=count where=(agree_eff=1) rename=(percent=pct_agree_eff));
run;


data &dsout. out.summary_se_&dsout;
 merge summary summary_disagree1 summary_eff_agree(drop=agree_eff);
 by true_pi;

 pct_agree_eff = pct_agree_eff / 100;

 /* if round(true_pi,1e-5) in (0.10,0.20,0.30,0.40,0.50) then do; */
  nDisp = strip(put(nObs,3.));
  nFinDisp  = strip(put(nFin,3.));
 /* end; */

 int_inc = 1 - (int_eff+int_fut);

 /*if round(true_pi,1e-5) in (0.10,0.20,0.30,0.40,0.50) then */ pDisp1 = put(int_eff,6.2);

 SS = strip(put(nObs,5.1))||' +  '||strip(put(nMis,5.1))||' = '||strip(put(nFin,5.1));
 MN = '(I) '||strip(put(int_mean,6.3))||' (F) '||strip(put(fin_mean,6.3));
 CP = '(I) '||strip(put(int_cover,6.3))||' (F) '||strip(put(fin_cover,6.3));
run;


proc means data = sim(rename=(fin_postProbSkept=finpp)) noprint nway;
 class true_pi;
 where int_eff=1;
 var agree_eff;
 output out = summary_finA(drop=_:) mean = agree_eff;
run;

proc means data = sim(rename=(fin_postProbSkept=finpp)) noprint nway;
 class true_pi;
 where int_eff=1 and agree_eff=0;
 var finpp;
 output out = summary_finB median= min= p1= p10= p25= p50= p75= p90= p99= max= / autoname;
run;

data summary_fin out.summary_fin_se_&dsout;
 merge summary_finA summary_finB;
 by true_pi;
run;


%mend sim_se;
