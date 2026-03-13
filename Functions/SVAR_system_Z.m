% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

function [b_T,D_p,A0_p,AL_p,BL_p,V_T,eps] = SVAR_system_Z(Y,P,X,S,D,ZeroRestrictions)

%Estimate and store the parameters at the grid of quantiles
[~,K] = size(Y);

[~,C] = size(D);

if isempty(X) == 0
    M = size(X,2);
else
    M = 0;
end

%Pre-allocate memory to store coefficients
b_T=[];
D_p = zeros(K,C);
A0_p = zeros(K,K);
if P > 0
    AL_p = zeros(K,K*P);
else
    AL_p = [];
end
if ismissing(X) == 0
    BL_p = zeros(K,M*(S+1));
else
    BL_p = [];
end
eps = [];

%Estimate system equation by equation
for i=1:K
    
    uq = [];
    
    %Estimation
    [b,~,u,~]= VAR_Estim_Z(Y,P,X,S,D,i,ZeroRestrictions);

    %Vector of all coefficients [(C+(j-1)*j/2+K*P+K*M*(S+1))*n*K x 1]
    b_T=[b_T;b];

    %Deterministic terms
    D_p(i,:) = b(1:C);

    %Simultaneous effects
    if i > 1
        A0_p(i,1:i-1) = b(end-(i-2):end);
    end

    if P > 0
        %Endogenous lags
        AL_p(i,:) = b(1+C:(C+K*P));
    end
    
    %Exogenous variables
    if isempty(X) == 0
        BL_p(i,:) = b(1+C+K*P:C+K*P+M*(S+1));
    end

    uq = [uq u];
    
    eps = [eps;uq];
    
end

%Cross-equation and -quantile covariance matrix
V_T = [];
