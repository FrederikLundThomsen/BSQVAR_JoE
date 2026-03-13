 ---------------------------------------------------------------- 
 ---------------------------------------------------------------- 
 Documentation for the BSQVAR Toolbox used in Chavleishvili, S., Engle, R., 
 Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
 Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
 A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.

 Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
 ---------------------------------------------------------------- 
 ---------------------------------------------------------------- 

 ---------------------------------------------------------------- 
 I) INTRODUCTION
 ---------------------------------------------------------------- 
 The BSQVAR toolbox replicates the analysis and figures from the above reference
 and is comprised of four folders:

 - Data

 - Figures

 - Functions

 - Scripts

 Each of these and their contents are described below.

 All codes can be run in Matlab, which also produces all relevant figures.

 ---------------------------------------------------------------- 
 II) DATA
 ---------------------------------------------------------------- 

 This folder contains two .csv-files: 
 i) US_data.csv
 ii) EA_data.csv

 The two files contain the data used to estimate the BSQVAR for the U.S.
 and euro area, respectively.

 They are directly loaded into Matlab using the relevant script

 ---------------------------------------------------------------- 
 III) FIGURES
 ---------------------------------------------------------------- 

 This folder contains replicates of all figures from the above reference.

 Each file is named after the figure number in the text.

 Running the relevant Matlab code will override the preexisting files
 unless otherwise specified.

 ---------------------------------------------------------------- 
 IV) FUNCTIONS
 ---------------------------------------------------------------- 

 This folder contains a series of Matlab support functions used to 
 replicate the analysis.

 They are:

 i)	QVAR_system_ZB.m
	- Draws from the posterior of the Bayesian SQVAR with potential zero restrictions

 ii)	igrnd.m
 	- Generate a draw from the inverse Gaussian distribution

 iii)	QVAR_Forecast_Scenario_Z_FixP.m
	- Projects the SQVAR forward taking the model parameters as input

 iv)	QVAR_Estim_Z.m
	- Frequentist univariate quantile regression at a single quantile with potential zero restrictions

 v)	QVAR_system_Z.m
	- Frequentist estimation of the SQVAR equation-by-equation for all quantiles with potential zero restrictions

 vi) 	bound.m 	
	- Part of the package of functions to perform frequentist quantile regression 

 vii) 	lp_fnm 		
	- Part of the package of functions to perform frequentist quantile regression 

 viii) 	rq_fnm.m 	
	- Part of the package of functions to perform frequentist quantile regression 

 ix)	SVAR_system_Z.m
	- Frequentist estimation of an SVAR with potential zero restrictions

 x)	VAR_system_Z.m
	- Frequentist VAR estimation with potential zero restrictions

 xi)	PlotFmt.m
	- Function to format an active Matlab figure

 ---------------------------------------------------------------- 
 V) SCRIPTS
 ---------------------------------------------------------------- 

 This folder contains the Matlab scripts used to replicate the analysis.

 Each script can, with some restrictions, be run individually or in complete sequence.

 The script 'S0_Wrapper.m' should be examined first for more information about how to run
 the different scripts as well as their outputs.

 ---------------------------------------------------------------- 
 ---------------------------------------------------------------- 
