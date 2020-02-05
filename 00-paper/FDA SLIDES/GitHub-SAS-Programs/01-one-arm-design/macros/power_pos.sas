

%macro plotA;

ods noproctitle;
title;
options topmargin=0.1in bottommargin=0.1in rightmargin=0.1in leftmargin=0.1in nodate nonumber;
option papersize=("9.7in","5.2in") nodate nonumber;

ods pdf file = "&outPath.\one-arm-design-characteristics-pos-A-&dsOut..pdf";
ods graphics  / height=5.0in width=9.5in noborder;
proc sgplot data = &dsout.;

 refline 0.20 0.40  / axis=x lineattrs=(thickness=2 pattern=2 color=verylightGray);

 yaxis values = (0 to 1 by 0.10) grid label='Probability';
 xaxis offsetmin=0.0625 offsetmax=0.0625 display=(nolabel);

 series x=true_pi y=int_eff / markers markerattrs=(symbol=circleFilled color=blue)  lineattrs=(color=blue pattern=1 thickness=2) datalabel=pDisp1;
 series x=true_pi y=int_fut / markers markerattrs=(symbol=squareFilled color=red)   lineattrs=(color=red  pattern=1 thickness=2);
 series x=true_pi y=int_inc / markers markerattrs=(symbol=diamondFilled color=gray) lineattrs=(color=gray pattern=1 thickness=2);

 xaxistable SS    / location=outside valueattrs=(Family=Arial Size=8 );
 xaxistable MN    / location=outside valueattrs=(Family=Arial Size=8 );
 xaxistable CP    / location=outside valueattrs=(Family=Arial Size=8 );



 label SS='SS' true_pi = 'True Response Probability';
 label int_eff = 'Stop Early for Efficacy'
       int_fut = 'Stop Early for Futility'
       int_inc = 'Inconclusive Findings w/ Full Dataset'
       MN = 'PM';

	   keylegend / position=top;
run;
ods pdf close;
%mend plotA;


%macro plotB;
ods noproctitle;
title;
options topmargin=0.1in bottommargin=0.1in rightmargin=0.1in leftmargin=0.1in nodate nonumber;
option papersize=("8.7in","5.2in") nodate nonumber;

ods pdf file = "&outPath.\one-arm-design-characteristics-pos-B-&dsOut..pdf";
ods graphics  / height=5.0in width=8.5in noborder;
ods graphics / attrpriority=none ;

proc sgplot data = summary_fin;
 band x=true_pi lower=finpp_min upper=finpp_max / fillattrs=(color=lightblue transparency=0.60) legendLabel='MIN-MAX' outline fill lineattrs=(color=lightblue);
 band x=true_pi lower=finpp_p1  upper=finpp_p99 / fillattrs=(color=blue transparency=0.80) legendLabel='P1-P99';
 
 band x=true_pi lower=finpp_p10 upper=finpp_p90 / fillattrs=(color=blue transparency=0.60) legendLabel='P10-P90';
 band x=true_pi lower=finpp_p25 upper=finpp_p75 / fillattrs=(color=blue transparency=0.40) legendLabel='P25-P75';
 series x=true_pi y=finpp_median / markers markerattrs=(symbol=diamondfilled color=veryDarkBlue) lineattrs=(color=veryDarkBlue) legendLabel='Median' datalabel=finpp_median;
 
 format agree_eff 6.2 finpp_median 6.2;
 xaxistable agree_eff / label='% AGR.' location=inside;
 refline 0.5 / axis=y lineattrs=(pattern=2 color=black);

 yaxis label = "Final P((*ESC*){unicode pi}>&PiNull.)" values=(0.50 to 1.0 by 0.10);
 xaxis label = 'True Response Rate' offsetmin=0.0625 offsetmax=0.0625;

 keylegend / position=top;
 
run;

ods pdf close;

%mend plotB;