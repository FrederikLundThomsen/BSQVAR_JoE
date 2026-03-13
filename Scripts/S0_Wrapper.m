% ---------------------------------------------------------------- %
% --------------------------- 0. Wrapper ------------------------- %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script replicates the analysis and figures using the BSQVAR as given
% in the above reference.

% All relevant scripts for different sections of the paper can be run from
% this unified wrapper.

% Supporting functions, estimation data and figures are provided within 
% this toolbox.

% NOTE: Scripts 1-3 are a prerequisite for the analysis and should 
%  always be run before any other script and in sequence.

% For further information, see also the README.txt file in the top level
% folder of this toolbox.

% ---------------------------------------------------------------- %

%Set working directory to location of scripts
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

% ---------------------------------------------------------------- %


%% ------------------------ PREREQUISITES ------------------------ %%

% I) Initialize the Matlab workspace
% ---------------------------------------------------------------- %
%    This script defines a series of global variables used througout.

S1_Initialization;

% ---------------------------------------------------------------- %

% II) Data
% ---------------------------------------------------------------- %
%    This script loads and transforms the data used for the analysis.

S2_Data;

% ---------------------------------------------------------------- %
       
% III) Estimation
% ---------------------------------------------------------------- %
%    This script estimates the BSQVAR equation-by-equation.
%
%    In addition, it produces the following figures from the paper:
%     - C.1: Euro area and U.S. time series
%     - D.1: Posterior inference for the euro area, omega(gamma)
%     - D.2: Posterior inference for the euro area, A0(gamma)
%     - D.3-D.6: Posterior inference for the euro area, A1(gamma)-A4(gamma)
%     - D.7: Posterior inference for the euro area exogenous variables, C(gamma)
%     - D.8: Posterior inference for the euro area dummy variables, B(gamma)
%     - D.9: Posterior means for the euro area tightness parameter, lambda(gamma)
%     - F.1: Posterior inference for the U.S., omega(gamma)
%     - F.2: Posterior inference for the U.S., A0(gamma)
%     - F.3-F.6: Posterior inference for the U.S., A1(gamma)-A4(gamma)
%     - F.7: Posterior inference for the U.S. exogenous variables, C(gamma)
%     - F.8: Posterior inference for the U.S. dummy variables, B(gamma)
%     - F.9: Posterior means for the U.S. tightness parameter, lambda(gamma)

S3_Estimation;

% ---------------------------------------------------------------- %

%% -------------------------- ANALYSIS -------------------------- %%

% IV) Produce posterior inference for quantile impulse response functions 
% (QIRFs).
% ---------------------------------------------------------------- %
%    This script obtains QIRFs through simulation methods as well as
%    posterior inference. See also 
%    - Chavleishvili, S. and Manganelli, S. (2023). 'Forecasting and 
%      Stress Testing with Quantile Vector Autoregressions'. 
%      Journal of Applied Econometrics, forthcoming 
%    - Chavleishvili, S., Kremer, M. and Lund-Thomsen, F. (2023).
%    'Quantifying financial stability trade-offs for monetary policy: a
%    quantile VAR approach'. ECB Working Paper Series, no. 2833.
%
%    In addition, it produces the following figures from the paper:
%     - 1: Quantile impulse response estimates
%     - F.10: U.S. quantile impulse response estimates

S4_QIRF;

% ---------------------------------------------------------------- %

% VI) Use the BSQVAR to compute counterfactual outcomes of different
% macroprudential policies.
% ---------------------------------------------------------------- %
%
%    In addition, it produces the following figures from the paper:
%     - 2: The financial cycle in two different experiments
%     - D.10: Posterior distributions of moments in counterfactual policy
%     scenario

S5_Counterfactual;

% ---------------------------------------------------------------- %

% VII) Use the BSQVAR to produce vintages of a metric for the 
% macroprudential policy stance under a passive and active policy scenario
% as well as posterior inference.
% ---------------------------------------------------------------- %
%
% NOTE: Executing the entire script is very time consuming and should only
% be done when necessary.
%
%    In addition, it produces the following figures from the paper:
%     - 3: Macro-prudential policy stance

S6_MacroprudentialPolicyStance;

%% END
