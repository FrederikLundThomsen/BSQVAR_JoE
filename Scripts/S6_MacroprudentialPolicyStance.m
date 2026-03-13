% ---------------------------------------------------------------- %
% ---------------- 6. Macroprudential policy stance -------------- %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script produces posterior inference for the macro-prudential policy
% stance and relates it to the macro-financial environment

% Figures: 
% - 3

% ---------------------------------------------------------------- %

%% ---- Macro-prudential policy stance

%Set forecasting horizon (H) and number of simulations (N) for each
%posterior draw
H = 4*4;
N = 50e+3;

%Initialize design matrixes
ImpulseDesign = zeros(H,K+M,nQ);
QuantileDesign = nan(H,K+M);
PathDesign = nan(H,K+M,nQ);

%Construct the baseline and policy scenarios

%Baseline is unconditional 
Baseline = {ImpulseDesign,QuantileDesign,PathDesign};

%Small shock to the financial cycle

%First obtain residuals at the median
[~,~,~,~,~,~,epsY] = QVAR_system_Z(Y_EA,P,X_EA,S,[D_EA(:,1) DM_ea],0.5,[ZeroRestrictions(1,:);ones(size(DM_ea,2),K);ZeroRestrictions(1+size(D_ea,2):end,:)]);
eps_std = std(reshape(epsY,[],K));

%Size of shock corresponding to a back-of-the-envelope estimate based on Ampudia, M., M. L. Duca, M. Farkas, G. Perez-Quiros, M. Pirovano, G. Runstler, and E. Tereanu (2021).
% On the effectiveness of macroprudential policy. ECB Discussion Paper 2559, 1-42
sri_shk = eps_std(sri_o)/(1*2);

%Shock period
shk_p = 1;
%Reversal period
rev_p = shk_p+ceil(H/2);

%Impose the shocks
ImpulseDesign(shk_p,sri_o,:) = -sri_shk;
ImpulseDesign(rev_p,sri_o,:) = sri_shk;

Policy = {ImpulseDesign,QuantileDesign,PathDesign};

%Number of draws from the posterior
pDraw = 2e+2;

%Threshold for growth shortfall/longrise
taug = 0;

%Index over vintages
PinT = P+1:1:Te;

%Arrays to store output
g_avg_t = zeros(length(PinT),H,2,pDraw);
g_shrt_t = zeros(length(PinT),H,2,pDraw);

%Start a timer
TIME = tic;

%Loop over the posterior draws for the euro area
for bp = 1:pDraw

    %Print progress
    "[" + string(datetime()) + "] " + "Time elapsed: " + sprintf("%s",duration([0, 0, toc(TIME)])) + " Iteration: " + bp + "."
    
    %Randomly select a posterior draw above the burn-in threshold
    rng(rnd_seed + bp)
    pstrr = randi([burn+1 biter],1,1);
    
    %Select the parameters given the posterior draw
    D_bb = D_ea(:,:,:,pstrr);
    A0_bb = A0_ea(:,:,:,pstrr);
    AL_bb = AL_ea(:,:,:,pstrr);
    if ~isempty(X_EA)
        BL_bb = BL_ea(:,:,:,pstrr);
        b_bb = b_bbx_EA(:,pstrr);
    else
        BL_bb = [];
        b_bb = [];
    end

    %For each vintage
    for t = 1:length(PinT)

        %Loop over the two scenarios
        for scen = 1:2

            if scen == 1
                ScenEval = Baseline;
            else
                ScenEval = Policy;
            end

            %Simulate the scenario with fiX_EAd parameters drawn from the posterior
            [Y_S,~] = QVAR_Forecast_Scenario_Z_FixP(Y_EA,P,X_EA,P,D_EA,Quant,PinT(t),H,N,ScenEval,ZeroRestrictions,D_bb,A0_bb,AL_bb,BL_bb,b_bb,rnd_seed);
        
            %Compute average growth and growth shortfall at each forecast
            %horizon

            %HxN matrix with all output paths
            g_cell = cell2mat(cellfun(@(c) c(2:end,gdp_o),Y_S,'UniformOutput',false)');

            %Take average over each horizon
            g_avg = mean(g_cell,2);
            g_avg_t(t,:,scen,bp) = g_avg;

            %Compute growth shortfall as the mean over all simulations,
            %where all values above the threshold, taug, are set equal to
            %zero
            g_shrt = mean(g_cell.*(g_cell < taug),2);
            g_shrt_t(t,:,scen,bp) = g_shrt;

        end

    end
    
end

%Parameters of policy function
%Discount factor
bdisc = 1;
%Weight on downside risks to growth
lambd = 1.5*2;
lambd_alt = 1.5*1;

%Quantify the policy maker's objective function for each policy
%scenario, at each vintage and for each posterior draw
ObjFun = zeros(length(PinT),2,pDraw);
ObjFunA = zeros(length(PinT),2,pDraw);
for bp = 1:pDraw
    for scen = 1:size(g_avg_t,3)
		ObjFun(:,scen,bp) = mean(bdisc.^(1:H).*(squeeze(g_avg_t(:,1:H,scen,bp)) + (lambd-1)*squeeze(g_shrt_t(:,1:H,scen,bp))),2);
		ObjFunA(:,scen,bp) = mean(bdisc.^(1:H).*(squeeze(g_avg_t(:,1:H,scen,bp)) + (lambd_alt-1)*squeeze(g_shrt_t(:,1:H,scen,bp))),2);
    end
end

%Take differences between the active(2) and passive(1) scenarios
ObjFunDiff = ObjFun(:,2,:) - ObjFun(:,1,:);
ObjFunDiffA = ObjFunA(:,2,:) - ObjFunA(:,1,:);

%% ---- Macro-financial environment

%Set up least squares regression with \delta u_{t} as the LHS 
% on macro-financial variables

%Number of RHS lags in addition to contemporaneous values
PP = 0;

%Compute the posterior mean of the macro prudential policy stance
uY = mean(ObjFunDiff,[2 3]);
%Define RHS
uX = rmmissing(lagmatrix([Y_EA X_EA],0));

%Limit the sample
uX = uX(P+1:end,:);

[EstCoeffCov,se,coeff] = hac(uX,100*uY,'Intercept',true,'Type','HAC','Weights','BT');

t_stat = coeff./se;

%% ---- Figures
%% ---- Figure 3: Macro-prudential policy stance

%HPDRs to plot
hpdr_out = 5;
hpdr_in = 32;

%Call a new figure
figure
tiledlayout(1,1,"TileSpacing","compact"',"Padding","compact")

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.25;
figp(4) = 21*3/5*1.125*1;
set(gcf,'Position',figp);

ylm = nan(2,1);
ylm(1) = min(ylm(1),min(prctile(ObjFunDiff,hpdr_out/2,[2 3])));
ylm(2) = min(ylm(2),max(prctile(ObjFunDiff,100-hpdr_out/2,[2 3])));
ylm = max(abs(ylm));

%Plot the posterior distribution of \delta u
nexttile
hold on
plot(qe_dt(PinT),mean(ObjFunDiff,[2 3]),'-','LineWidth',2,'Color',dflt_col{1})
plot(qe_dt(PinT(1)),nan,'-.','LineWidth',2,'Color',dflt_col{2})
patch([(qe_dt(PinT))' flip(qe_dt(PinT))'],[prctile(ObjFunDiff,hpdr_in/2,[2 3])' flip(prctile(ObjFunDiff,100-hpdr_in/2,[2 3]))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.4,'LineWidth',1)
patch([(qe_dt(PinT))' flip(qe_dt(PinT))'],[prctile(ObjFunDiff,hpdr_out/2,[2 3])' flip(prctile(ObjFunDiff,100-hpdr_out/2,[2 3]))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1)
plot(qe_dt(PinT),mean(ObjFunDiffA,[2 3]),'-.','LineWidth',2,'Color',dflt_col{2})
yline(0)
hold off
axis tight
ylim([-0.04 0.04])
box on
ylabel('Too tight / Too loose','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
title('$\lambda^{p} = 1.5$','interpreter','latex')

%Define the legend
lh = legend([ ...
      "$\Delta u$ for $\lambda^{p}=" + string(sprintf('%1.1f',lambd)) + "$ (Posterior mean with " + (100 - hpdr_in) + "\% and " + (100 - hpdr_out) + "\% credible bands)"...
      ,"$\Delta u$ for $\lambda^{p}=" + string(sprintf('%1.1f',lambd_alt)) + "$ (Posterior mean)"...
    ],'Orientation','Horizontal','Location','Northoutside','NumColumns',4,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Add euro area recessions as defined by CEPR
rc = qe_dt >= datetime('1992-03-31','Format','yyyy-MM-dd') & qe_dt < datetime('1993-10-01','Format','yyyy-MM-dd');
patch([qe_dt(rc)' flip(qe_dt(rc))'],[-15*ones(1,sum(rc)),15*ones(1,sum(rc))],dflt_col{16},'FaceAlpha',0.25,'EdgeColor','none')
hold on
rc = qe_dt >= datetime('2008-03-31','Format','yyyy-MM-dd') & qe_dt < datetime('2009-07-01','Format','yyyy-MM-dd');
patch([qe_dt(rc)' flip(qe_dt(rc))'],[-15*ones(1,sum(rc)),15*ones(1,sum(rc))],dflt_col{16},'FaceAlpha',0.25,'EdgeColor','none')
hold on
rc = qe_dt >= datetime('2011-09-30','Format','yyyy-MM-dd') & qe_dt < datetime('2013-04-01','Format','yyyy-MM-dd');
patch([qe_dt(rc)' flip(qe_dt(rc))'],[-15*ones(1,sum(rc)),15*ones(1,sum(rc))],dflt_col{16},'FaceAlpha',0.25,'EdgeColor','none')
hold on
rc = qe_dt >= datetime('2019-12-31','Format','yyyy-MM-dd') & qe_dt < datetime('2020-07-01','Format','yyyy-MM-dd');
patch([qe_dt(rc)' flip(qe_dt(rc))'],[-15*ones(1,sum(rc)),15*ones(1,sum(rc))],dflt_col{16},'FaceAlpha',0.25,'EdgeColor','none')
hold off

%Ad hoc modifications
fz = gcf();
fz.Children(1).Children(1).String = fz.Children(1).Children(1).String(1:2);
fz.Children(1).Children(1).FontSize = fz.Children(1).Children(1).FontSize + 3;
fz.Children(1).Children(2).FontSize = fz.Children(1).Children(2).FontSize + 2;

fz.Children(1).Children(2).YAxis(2).TickValues = [];

%Print
print(fig_pth+"\MAIN_3",'-dpdf','-vector','-bestfit');

%% END