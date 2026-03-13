% ---------------------------------------------------------------- %
% ---------------------------- 2. DATA --------------------------- %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script loads and transforms the data used in the analysis

% Figures: ~

% ---------------------------------------------------------------- %

%% ---- US data

%Load the data
US_dt_q = readtimetable(dt_pth + "\US_data.csv","RowTimes",'Time');

%Transform the data as needed
US_dt_q.gdp_d = 400*(log(US_dt_q.REALGDP) - log(lagmatrix(US_dt_q.REALGDP,1)));
US_dt_q.cpi_d = 400*(log(US_dt_q.CPI) - log(lagmatrix(US_dt_q.CPI,1)));
US_dt_q.spg_d = 400*(log(US_dt_q.SPGSCISADJ) - log(lagmatrix(US_dt_q.SPGSCISADJ,1)));
US_dt_q.int_l = US_dt_q.FEDFUNDS;

%Collect estimation sample
Yu = [US_dt_q.CISS US_dt_q.SRI US_dt_q.cpi_d US_dt_q.gdp_d US_dt_q.int_l];

%Exogenous variables
X_US = [US_dt_q.spg_d];

%Vector with time periods
US_q = US_dt_q.Time;

%Remove observations with NaNs
if isempty(X_US)
    nan_rm = any(isnan(Yu),2);
else
    nan_rm = any(isnan([Yu X_US]),2);
end
US_qe = US_q(~nan_rm);
Yu = Yu(~nan_rm,:);
if ~isempty(X_US)
    X_US = X_US(~nan_rm,:);
end

%Construct deterministic terms
D_US = ones(size(Yu,1),1);
%   Dummies for COVID19
DM_US = [];
DM_US = [DM_US double(logical(year(US_qe) == 2020 & quarter(US_qe) == 1))];
DM_US = [DM_US double(logical(year(US_qe) == 2020 & quarter(US_qe) == 2))];
DM_US = [DM_US double(logical(year(US_qe) == 2020 & quarter(US_qe) == 3))];
DM_US = [DM_US double(logical(year(US_qe) == 2020 & quarter(US_qe) == 4))];

%Gather them in one matrix
D_US = [D_US DM_US];

%% ---- Euro area data

%Load the data
EA_dt_q = readtimetable(dt_pth + "\EA_data.csv","RowTimes",'Time');

%Transform the data as needed
EA_dt_q.gdp_d = 400*(log(EA_dt_q.REALGDP) - log(lagmatrix(EA_dt_q.REALGDP,1)));
EA_dt_q.cpi_d = 400*(log(EA_dt_q.HICP) - log(lagmatrix(EA_dt_q.HICP,1)));
EA_dt_q.spg_d = 400*(log(EA_dt_q.SPGSCISADJ) - log(lagmatrix(EA_dt_q.SPGSCISADJ,1)));
EA_dt_q.int_l = EA_dt_q.OIS3M;

%Collect estimation sample
Ye = [EA_dt_q.CISS EA_dt_q.SRI EA_dt_q.cpi_d EA_dt_q.gdp_d EA_dt_q.int_l];

%Exogenous variables
X_EA = [EA_dt_q.spg_d];

%Vector with time series periods
q_dt = EA_dt_q.Time;

%Remove observations with NaNs
if isempty(X_EA)
    nan_rm = any(isnan(Ye),2);
else
    nan_rm = any(isnan([Ye X_EA]),2);
end
qe_dt = q_dt(~nan_rm);
Ye = Ye(~nan_rm,:);
if ~isempty(X_EA)
    X_EA = X_EA(~nan_rm,:);
end

%Deterministic terms
D_EA = ones(size(Ye,1),1);
%   Dummies for COVID19
DM_ea = [];
DM_ea = [DM_ea double(logical(year(qe_dt) == 2020 & quarter(qe_dt) == 1))];
DM_ea = [DM_ea double(logical(year(qe_dt) == 2020 & quarter(qe_dt) == 2))];
DM_ea = [DM_ea double(logical(year(qe_dt) == 2020 & quarter(qe_dt) == 3))];
DM_ea = [DM_ea double(logical(year(qe_dt) == 2020 & quarter(qe_dt) == 4))];

D_EA = [D_EA DM_ea];

%Data dimensions
C = size(D_EA,2);
K = size(Ye,2);
M = size(X_EA,2);

%% END