
/*************************************************************************/
/**************         Macros for Reference interval          **************/
/************************     Author:Xiao Ma        **********************/
/*************************************************************************/
    
/*************************************************************************/
/**************         Data Imputing and Cleaning          **************/
/*************************************************************************/
libname PlanB "/home/u59750830/PLAN B";

proc import datafile="/home/u59750830/PLAN B/healthy2.csv" dbms=csv 
		out=healthy2 replace;
run;

/* creat new variables */
DATA healthy2;
	SET healthy2;
	log_lsm=log(stiff_measure);
	StudyID=s;
RUN;

/* 95% confidence intervals for the study means */
proc means data=healthy2 std mean noprint;
	var stiff_measure;
	class StudyID;
	ways 1;
	output out=summary_data std=SD mean=mean;
run;

/* Deriving Aggregate Dataset */
data PLANB.Aggregate_data;
	set SUMMARY_DATA;
	lower_mean=mean-quantile("t", 0.975, _FREQ_-1)* SD/SQRT(_FREQ_);
	upper_mean=mean+quantile("t", 0.975, _FREQ_-1)* SD/SQRT(_FREQ_);
	n=_FREQ_;
	mean_bak=mean;
	logmean=log(mean/sqrt(1 + ((n-1)/(n))*sd**2/mean**2));
	logsd=sqrt(log(1 + ((n-1)/(n))*sd**2/mean**2));
	keep StudyID sd n lower_mean mean mean_bak upper_mean lower_mean logmean logsd;
run;

proc export data=PLANB.Aggregate_data 
		outfile="/home/u59750830/PLAN B/aggregate_data.csv" dbms=csv replace;
run;

/*************************************************************************/
/******************        Define Global Macro          ******************/
/*************************************************************************/

/* Normal_calculation */
%macro normal_calculation;
	data result;
		lower=quantile('NORMAL', 0.025, &mean, &sd);
		mean=&mean;
		upper=quantile('NORMAL', 0.975, &mean, &sd);
	run;

%mend;

/* Lognormal_calculation */
%macro lognormal_calculation;
	data result;
		lower=exp(quantile('NORMAL', 0.025, &mean, &sd));
		mean=exp(&mean);
		upper=exp(quantile('NORMAL', 0.975, &mean, &sd));
	run;

%mend;

/*************************************************************************/
/*****************************   METHOD 1  *******************************/
/*                      frequentist approach for AD                      */
/*************************************************************************/
%macro FreqFit(data_use=, mean_use=, id_use=, sd_use=, n_use=, dist=);
	data use;
		set &data_use;
		var=&sd_use**2/&n_use;
		col=_n_;
		row=_n_;
		value=var;
	run;

	ods trace on;
	ods output CovParms=CovParms_AD_freq;
	ods output SolutionF=SolutionF_AD_freq;
	ods exclude ModelInfo Dimensions NObs CovParms FitStatistics SolutionF;

	proc mixed data=use order=data method=reml NOCLPRINT NOITPRINT;
		class &id_use;
		model &mean_use= /solution cl ddfm=kr;
		random &id_use/gdata=use;
		repeated diag;
	run;

	ods trace off;

	/* create macro tau_hat */
	data _null_;
		set CovParms_AD_freq;
		tau_hat=sqrt(Estimate);
		call symput("tau_hat", tau_hat);
	run;

	%put tau_hat = &tau_hat;

	/* create macro sigma_hat */
	proc iml;
		varNames={"&mean_use" "&sd_use" "&n_use"};
		use &data_use;
		read all var varNames;
		sigma_hat=sqrt((&n_use-1)`*(&sd_use##2)/sum(&n_use-1));
		call symputx("sigma_hat", rowcat(char(sigma_hat)));
		%put &sigma_hat;
		quit;

	/* create macro pooled_mean */
	data _null_;
		set SolutionF_AD_freq;
		pooled_mean=Estimate;
		sd=sqrt(&sigma_hat**2+&tau_hat**2);
		call symput("mean", pooled_mean);
		call symput("sd", sd);
	run;

	%put mean = &mean;
	%put sd = &sd;

	/*log normal distributed */	
	%if &dist=lognormal %then
		%do;
			%lognormal_calculation;
		%end;

	/* normal distributed */
    %else
		%do;
			%normal_calculation;
		%end;

	proc print data=result NOOBS;
		title1 "Reference interval for AD using frequentist method";
		title2 "&data_use, &mean_use with &dist distribution";
	run;

%mend;

/*************************************************************************/
/*****************************   METHOD 2  *******************************/
/*                      frequentist approach for IPD                     */
/*************************************************************************/
%macro FreqIPD(data_use=, variable_use=, id_use=, sd_use=, n_use=, dist=);

	ods trace on;
	ods output SolutionF=SolutionF_IPD_freq;
	ods output CovParms=CovParms_IPD_freq;	
	ods exclude ModelInfo ClassLevels Dimensions NObs IterHistory 
		ConvergenceStatus CovParms FitStatistics SolutionF;
	proc mixed data=&data_use;
		class &id_use;
		model &variable_use= / s cl;
		random intercept / subject=&id_use;
	run;

	ods trace off;

	/* pooled mean */
	data _null_;
		set SolutionF_IPD_freq;
		mean=Estimate;
		call symput("mean", mean);
	run;

	%put mean = &mean;

	/* sigma and tau */
	data _null_;
		set CovParms_IPD_freq;
		sigma_IPD_freq=Estimate;
		call symput("sigma_IPD_freq", sigma_IPD_freq);
	run;

	%put sigma_IPD_freq = &sigma_IPD_freq;

	data _null_;
		set CovParms_IPD_freq;
		where subject="&id_use";
		sd=sqrt(Estimate+&sigma_IPD_freq);
		call symput("sd", sd);
	run;

	%put sd = &sd;

	/*log normal distributed */
	%if &dist=lognormal %then
		%do;
			%lognormal_calculation;
		%end;

	/* normal distributed */
    %else
		%do;
			%normal_calculation;
		%end;

	proc print data=result NOOBS;
		title1 "Reference interval for IPD using frequentist method";
		title2 "&data_use, &variable_use with &dist distribution";
	run;

%mend;

/*************************************************************************/
/*****************************   METHOD 3  *******************************/
/*                        Empirical approach for AD                      */
/*************************************************************************/
%macro EmpFit(data_use=, mean_use=, id_use=, sd_use=, n_use=, dist=);
	proc iml;
		varNames={"&mean_use" "&sd_use" "&n_use"};
		use &data_use;
		read all var varNames;
		sigma_hat=sqrt((&n_use-1)`*(&sd_use##2)/sum(&n_use-1));
		mean=&n_use`*&mean_use/sum(&n_use);
		var_ew=(&n_use-1)`*(&mean_use-mean)##2/sum(&n_use-1);
		sd=sqrt(var_ew+sigma_hat##2);
		call symputx("sigma_hat", rowcat(char(sigma_hat)));
		call symputx("mean", rowcat(char(mean)));
		call symputx("var_ew", rowcat(char(var_ew)));
		call symputx("sd", rowcat(char(sd)));
		%put &sigma_hat;
		%put &mean;
		%put &var_ew;
		%put &sd;
	quit;

		/*log normal distributed */
		%if &dist=lognormal %then
			%do;
				%lognormal_calculation;
			%end;

		/* normal distributed */
        %else
			%do;
				%normal_calculation;
			%end;

	proc print data=result NOOBS;
		title1 "Reference interval for AD using empirical method";
		title2 "&data_use, &mean_use with &dist distribution";
	run;

%mend;

/*************************************************************************/
/*****************************   METHOD 4  *******************************/
/*                        Empirical approach for IPD                     */
/*************************************************************************/
%macro EmpIPD(data_use=, variable_use=, dist=);
	proc iml;
		varNames={"&variable_use"};
		use &data_use;
		read all var varNames;

		%if &dist=lognormal %then
			%do;
				mean_ov=mean(log(&variable_use));
				sd_ov=std(log(&variable_use));
				lower=exp(quantile("normal", 0.025, mean_ov, sd_ov));
				upper=exp(quantile("normal", 0.975, mean_ov, sd_ov));
			%end;
		%else
			%do;
				lower=quantile("normal", 0.025, mean(&variable_use), std(&variable_use));
				upper=quantile("normal", 0.975, mean(&variable_use), std(&variable_use));
			%end;
		call symputx("lower", rowcat(char(lower)));
		call symputx("upper", rowcat(char(upper)));
		%put &lower &upper;

	data result;
		lower=&lower;
		upper=&upper;
	run;

	proc print data=result NOOBS;
		title1 "Reference interval for IPD using empirical method";
		title2 "&data_use, &variable_use with &dist distribution";
	run;

%mend;

/*************************************************************************/
/*****************************   METHOD 5  *******************************/
/*                  Mixture distribution approach for IPD                */
/*************************************************************************/
%macro weight_quantile;
	varNames={"&mean_use" "&sd_use" "&n_use"};
	use &data_use;
	read all var varNames;

	%if &weight=ss %then
		%do;
			weight=&n_use/sum(&n_use);
		%end;
	%else
		%do;
			weight=(&n_use/&sd_use##2)/sum(&n_use/&sd_use##2);
		%end;

	%if &dist=normal %then
		%do;
			quantile=cdf("normal", x, &mean_use, &sd_use);
		%end;
	%else %if &dist=lognormal %then
		%do;
			meanlog=log(&mean_use/sqrt(1+(&sd_use##2)/(&mean_use##2)));
			sdlog=sqrt(log(1+(&sd_use##2/&mean_use##2)));
			quantile=cdf("lognormal", x, meanlog, sdlog);
		%end;
	%else
		%do;
			shape=&mean_use##2/sd##2;
			scale=sd##2/&mean_use;
			quantile=cdf("gamma", x, shape, scale);
		%end;
	sum1=quantile`*weight;
	sum=sum1[1];
%mend;

%macro MixFit(data_use=, mean_use=, sd_use=, n_use=, lower_quantile=, 
		upper_quantile=, dist=, weight=);
	proc iml;
		/* Lower bound*/
		start Lower(x);
		%weight_quantile;
		return(sum-&lower_quantile);
		finish;
		intervals={-100 100};
		lower=froot("Lower", intervals);
		call symputx("lower", rowcat(char(lower)));
		%put &lower;

		/* Upper bound*/
		start Upper(x);
		%weight_quantile;
		return(sum-&upper_quantile);
		finish;
		Upper=froot("Upper", intervals);
		call symputx("upper", rowcat(char(upper)));
		%put &upper;

	data result;
		lower=&lower;
		upper=&upper;
	run;

	proc print data=result NOOBS;
		title1 "Reference interval for AD using mix distribution method";
		title2 "&data_use, &mean_use with &dist distribution";
		title3 "Reference interval [&lower_quantile ~ &upper_quantile], weight=&weight";
	run;

%mend;
	
/*************************************************************************/
/*************************    MACRO TEST      ****************************/
/*************************************************************************/
%FreqFit(data_use=PLANB.Aggregate_data, mean_use=logmean, id_use=studyID, sd_use=logsd, n_use=n, dist=lognormal) 

%FreqIPD(data_use=healthy2, variable_use=stiff_measure, id_use=s, dist=normal) 

%EmpFit(data_use=PLANB.Aggregate_data, mean_use=mean, id_use=studyID, sd_use=sd, n_use=n, dist=normal) 

 
	
%MixFit(data_use=PLANB.Aggregate_data, mean_use=mean, sd_use=SD, n_use=n, lower_quantile=0.025, upper_quantile=0.975, dist=normal, weight=ss) 

%MixFit(data_use=PLANB.Aggregate_data, mean_use=mean, sd_use=SD, n_use=n, lower_quantile=0.025, upper_quantile=0.975, dist=lognormal, weight=tt) 

%MixFit(data_use=PLANB.Aggregate_data, mean_use=mean, sd_use=SD, n_use=n, lower_quantile=0.025, upper_quantile=0.975, dist=gamma, weight=tt)
