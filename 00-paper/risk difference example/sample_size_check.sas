data check;
call streaminit(1);
  n = 112;
   pi0 = 0.40;
   pi1 = %sysevalf(0.4 + (0.67-0.40)*0.50);
   nSim = 1000000;
   do sim = 1 to nSim;
     y = rand('binomial',pi1,n);
     r + (sdf('beta',pi0,y+0.5,(n-y)+0.5)>0.975);
   end;
   r=r/nSim;
   put "power=" r;
run;

proc power;
   onesamplefreq test=exact
      nullproportion = 0.4
      proportion = 0.67
      ntotal = 1 to 100 by 1
      power = .;
run;

proc power;
   onesamplefreq test=exact
      nullproportion = 0.4
      proportion = 0.535
      ntotal = 1 to 100 by 1
      power = .;
run;

proc power;
   onesamplefreq test=z method=normal
      nullproportion = 0.4
      proportion = 0.535
      sides = 2
      ntotal = .
      power = .8;
run;
