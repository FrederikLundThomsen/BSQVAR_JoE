% ---------------------------------------------------------------- %
% ----------------------- 1. INITIALIZATION ---------------------- %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script initializes the Matlab environment

% Figures: ~

% ---------------------------------------------------------------- %

%Clear workspace and close all open windows
clear
clc
close all

%Point to data location
dt_pth = strrep(cd,'\Scripts','\Data');

%Path to dump figure prints
fig_pth = strrep(cd,'\Scripts','\Figures');

%Point to central Matlab functions
addpath(strrep(cd,'\Scripts','\Functions'));

%List of default color codes for plotting
dflt_col = {
    [0,50,153]/255,...
    [255,180,0]/255,...
    [255,75,0]/255,...
    [101,184,0]/255,...
    [0,177,234]/255,...
    [0,120,22]/255,...
    [129,57,198]/255,...
    [92,92,92]/255,...
    [152,161,208]/255,...
    [253,221,167]/255,...
    [246,177,131]/255,...
    [206,225,175]/255,...
    [215,238,248]/255,...
    [141,184,141]/255,...
    [174,151,199]/255,...
    [169,169,169]/255,...
    [217,217,217]/255
    };

%Seed for the random number generation
rnd_seed = 17071992;

%% END