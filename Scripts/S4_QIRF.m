% ---------------------------------------------------------------- %
% --------------------------- 4. QIRF ---------------------------- %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script produces the quantile impulse response functions for the euro
% area and US BSQVARs

% Figures: 
% - 1
% - H.10

% ---------------------------------------------------------------- %

%% ---- Euro area QIRFs

%Set forecasting horizon (H) and number of simulations (N)
H = 20;
N = 20e+3;

%Specify the quantiles of interest
Qplot = [find(round(Quant,2) == 0.1) find(round(Quant,2) == 0.5) find(round(Quant,2) == 0.9)];

%Harmonize forecast origin for the forward simulations
%We use the median values over the sample period
InitCndY = repmat(prctile(Y_EA,50),P+1,1);
InitCndD = repmat([1 zeros(1,C-1)],P+1,1);
if ~isempty(X_EA)
    InitCndX = repmat(prctile(X_EA,50),P+1,1);
else
    InitCndX = [];
end

%Specify shocks to the system, one shock for each Y and X
%We use the standard deviation of model residuals at the frequentist 
%median QVAR regression
[~,~,~,~,~,~,epsY] = QVAR_system_Z(Y_EA,P,X_EA,S,[D_EA(:,1) DM_ea],0.5,[ZeroRestrictions(1,:);ones(size(DM_ea,2),K);ZeroRestrictions(1+size(D_ea,2):end,:)]);
sigma = std(reshape(epsY,[],K));

if isempty(X_EA) == 0
    for m = 1:M
        [~,~,~,~,~,~,epsX] = QVAR_system_Z(X_EA,P,[],0,D_EA,0.5,[]);
        sigma = [sigma std(epsX)];
    end
end

%Number of draws from the posterior
pDraw = 4e+2;

%Initialize array to store output
QIRF_ea = zeros(K,K+M,H+1,size(Qplot,2),pDraw);

%Loop over the posterior draws for the euro area
for bp = 1:pDraw

    if bp == 1
        TIME = tic;
    end

    %Print progress
    "[" + string(datetime()) + "] " + "Time elapsed: " + sprintf("%s",duration([0, 0, toc(TIME)])) + " Iteration: " + bp + "."
    
    %Randomly select a posterior draw above the burn-in threshold
    rng(rnd_seed + bp)
    pstrr = randi([burn+1 biter],1,1);
    
    %Select the model parameters from the given draw
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

    %Simulate baseline path, cf. Chavleishvili et al. 2021 and 
    %Chavleishvili et al. 2023 with fixed parameters drawn from the 
    %posterior

    %Note that the estimation sample is stacked on top of the specified 
    %initial conditions and the forecast origin shifted forward accordingly
    [Y_F,~] = QVAR_Forecast_Scenario_Z_FixP([Y_EA;InitCndY],P,[X_EA;InitCndX],P,[D_EA;InitCndD],Quant,size(Y_EA,1)+P+1,H,N,{[],[],[]},ZeroRestrictions,D_bb,A0_bb,AL_bb,BL_bb,b_bb,rnd_seed);
    
    %Initialize temporary array to store the non-shocked model projection
    q_base = zeros(K,H+1,size(Qplot,2));
    for k=1:K
        q_base(k,:,:)=prctile(cell2mat(cellfun(@(c)c(:,k),Y_F,'UniformOutput',false)'),Quant(Qplot)*100,2);
    end
    
    %Loop over each shock
    for a=1:K+M
            
        %Construct design matrix for shocks to endogenous and exogenous variables
        Shocks = zeros(H,K+M,nQ);
    
        %Add the shock
        Shocks(1,a,:) = sigma(a);
    
        %Simulate paths for the given specification of shocks
        [Y_S,~] = QVAR_Forecast_Scenario_Z_FixP([Y_EA;InitCndY],P,[X_EA;InitCndX],P,[D_EA;InitCndD],Quant,size(Y_EA,1)+P+1,H,N,{Shocks,[],[]},ZeroRestrictions,D_bb,A0_bb,AL_bb,BL_bb,b_bb,rnd_seed);
    
        %Obtain the IRFs for the full set of simulations for each quantile
        %and store them
        q_shock = zeros(K,H+1,size(Qplot,2));
        for k=1:K
            q_shock(k,:,:)=prctile(cell2mat(cellfun(@(c)c(:,k),Y_S,'UniformOutput',false)'),Quant(Qplot)*100,2);
            
            %Compute the QIRF as the difference between the shocked and
            %non-shocked model projection
            QIRF_ea(k,a,:,:,bp) = squeeze(q_shock(k,:,:))-squeeze(q_base(k,:,:));
        end
        clearvars Y_S
    end

end

%% ---- US QIRFs

%Initialize array to store output
QIRF_us = zeros(K,K+M,H+1,size(Qplot,2),pDraw);

%Harmonize forecast origin for the forward simulations
%We use the median values over the sample period
InitCndY_US = repmat(prctile(Y_US,50),P,1);
InitCndD_US = repmat([1 zeros(1,C-1)],P,1);
if ~isempty(X_US)
    InitCndX_US = repmat(prctile(X_US,50),P,1);
else
    InitCndX_US = [];   
end

%Specify shocks to the system, one shock for each Y and X
%We use the standard deviation of model residuals at the frequentist 
%median QVAR regression
[~,~,~,~,~,~,epsY] = QVAR_system_Z(Y_US,P,X_US,S,[D_US(:,1) DM_US],0.5,[ZeroRestrictions(1,:);ones(size(DM_US,2),K);ZeroRestrictions(1+size(D_US,2):end,:)]);
sigmaUS = std(reshape(epsY,[],K));

if isempty(X_US) == 0
    for m = 1:M
        [~,~,~,~,~,~,epsX] = QVAR_system_Z(X_US,P,[],0,D_US,0.5,[]);
        sigmaUS = [sigmaUS std(epsX)];
    end
end

%Loop over the posterior draws for the US
for bp = 1:pDraw

    if bp == 1
        TIME = tic;
    end

    %Print progress
    "[" + string(datetime()) + "] " + "Time elapsed: " + sprintf("%s",duration([0, 0, toc(TIME)])) + " Iteration: " + bp + "."
    
    %Randomly select a posterior draw above the burn-in threshold
    rng(rnd_seed + bp)
    pstrr = randi([burn+1 biter],1,1);
    
    %Select the model parameters from the given draw
    D_bb = D_us(:,:,:,pstrr);
    A0_bb = A0_us(:,:,:,pstrr);
    AL_bb = AL_us(:,:,:,pstrr);
    if ~isempty(X_US)
        BL_bb = BL_us(:,:,:,pstrr);
        b_bb = b_bbx_US(:,pstrr);
    else
        BL_bb = [];
        b_bb = [];
    end

    %Simulate baseline path, cf. Chavleishvili et al. 2021 and 
    %Chavleishvili et al. 2023 with fixed parameters drawn from the 
    %posterior

    %Note that the estimation sample is stacked on top of the specified 
    %initial conditions and the forecast origin shifted forward accordingly
    [Y_F,~] = QVAR_Forecast_Scenario_Z_FixP([Y_US;InitCndY_US],P,[X_US;InitCndX_US],P,[D_US;InitCndD_US],Quant,size(Y_US,1)+P,H,N,{[],[],[]},ZeroRestrictions,D_bb,A0_bb,AL_bb,BL_bb,b_bb,rnd_seed);
    
    %Initialize temporary array to store the non-shocked model projection
    q_base = zeros(K,H+1,size(Qplot,2));
    for k=1:K
        q_base(k,:,:)=prctile(cell2mat(cellfun(@(c)c(:,k),Y_F,'UniformOutput',false)'),Quant(Qplot)*100,2);
    end
            
    %Loop over each shock
    for a=1:K+M
            
        %Construct design matrix for shocks to endogenous and exogenous variables
        Shocks = zeros(H,K+M,nQ);
    
        %Add the shock
        Shocks(1,a,:) = sigmaUS(a);
    
        %Simulate paths for the given specification of shocks
        [Y_S,~] = QVAR_Forecast_Scenario_Z_FixP([Y_US;InitCndY_US],P,[X_US;InitCndX_US],P,[D_US;InitCndD_US],Quant,size(Y_US,1)+P,H,N,{Shocks,[],[]},ZeroRestrictions,D_bb,A0_bb,AL_bb,BL_bb,b_bb,rnd_seed);
    
        %Obtain the IRFs for the full set of simulations for each quantile
        %and store them
        q_shock = zeros(K,H+1,size(Qplot,2));
        for k=1:K
            q_shock(k,:,:)=prctile(cell2mat(cellfun(@(c)c(:,k),Y_S,'UniformOutput',false)'),Quant(Qplot)*100,2);
            %Compute the QIRF as the difference between the shocked and
            %non-shocked model projection
            QIRF_us(k,a,:,:,bp) = squeeze(q_shock(k,:,:))-squeeze(q_base(k,:,:));
        end
        clearvars Y_S
    end
end

%% ---- Figures

alpha = 5;

%% ---- Figure 1: Quantile impulse response functions for the euro area

%Call a new figure
figure
tiledlayout(K,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 29.7;
figp(4) = 21*0.9;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','landscape')

%Loop over endogenous variables
for k=1:K
    %Loop over shocks
    for a=1:K
        nexttile
        ymx = 0;
        ymn = 0;
        for q = 1:length(Qplot)
            plot([0:H-1],squeeze(mean(QIRF_ea(k,a,2:end,q,:),[1 2 5])),'-','Color',dflt_col{q},'LineWidth',2)
            hold on
        end
        ymx = max([ymx;max(squeeze(mean(QIRF_ea(k,a,2:end,:,:),[1 2 5])),[],'all')]);
        ymn = min([ymn;min(squeeze(mean(QIRF_ea(k,a,2:end,:,:),[1 2 5])),[],'all')]);
        %Loop over quantiles of interest
        for q=1:size(Qplot,2)
            patch([(0:H-1) fliplr(0:H-1)],[squeeze(prctile(QIRF_ea(k,a,2:end,q,:),alpha/2,[1 2 4 5]))' fliplr(squeeze(prctile(QIRF_ea(k,a,2:end,q,:),100-alpha/2,[1 2 4 5]))')],dflt_col{q},'EdgeColor',dflt_col{q},'LineStyle','--','FaceAlpha',.2)
            ymx = max([ymx;max(squeeze(prctile(QIRF_ea(k,a,2:end,q,:),100-alpha/2,[1 2 4 5])))]);
            ymn = min([ymn;min(squeeze(prctile(QIRF_ea(k,a,2:end,q,:),alpha/2,[1 2 4 5])))]);
            hold on
        end
        yline(0)
        ylm = max(abs(ymx),abs(ymn));
        axis tight; box on;
        hold off
        ylim([-ylm ylm])
        xlim([0 H-1])
        xticks((0:4:H))
        ax = gca();
        Yzr_pos = find(ax.YTick == 0);
        if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
            ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
        end
        if k == K
            xlabel('Quarters ahead','interpreter','latex')
        end
        if k == 1
            if a == gdp_o
                title("Shock to real GDP growth",'interpreter','latex');
            elseif a == sri_o
                title("Shock to the fin. cycle",'interpreter','latex');                
            else
                title("Shock to "+ttl{a},'interpreter','latex');
            end
        end
        if a == 1
                if k == irt_o
                    ylabel('3M OIS rate','interpreter','latex');
                elseif k == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif k == inf_o
                    ylabel(["HICP";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{k}],'interpreter','latex');
                end
        end
        set(gca,'TickLabelInterpreter','latex')
    end
end

%Define the legend
lh = legend(...
    round(Quant(Qplot),2) + " quantile (Posterior mean + 95\% credible band)"...
    ,'Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for c = 2:length(fz.Children(1).Children)
    fz.Children.Children(c).YAxis(2).TickValues = [];

    fz.Children.Children(c).FontSize = fz.Children.Children(c).FontSize + 1;
end

fz.Children.Children(1).FontSize = fz.Children.Children(1).FontSize + 1;

print(fig_pth + "\MAIN_1",'-dpdf','-vector','-bestfit');

%% ---- ---- Web appendix
%% ---- ---- Figure H.10: U.S. quantile impulse response functions

%Call a new figure
figure
tiledlayout(K,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 29.7;
figp(4) = 21*0.9;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','landscape')

%Loop over endogenous variables
for k=1:K
    %Loop over shocks
    for a=1:K
        nexttile
        ymx = 0;
        ymn = 0;
        for q = 1:length(Qplot)
            plot([0:H-1],squeeze(mean(QIRF_us(k,a,2:end,q,:),[1 2 5])),'-','Color',dflt_col{q},'LineWidth',2)
            hold on
        end
        ymx = max([ymx;max(squeeze(mean(QIRF_us(k,a,2:end,:,:),[1 2 5])),[],'all')]);
        ymn = min([ymn;min(squeeze(mean(QIRF_us(k,a,2:end,:,:),[1 2 5])),[],'all')]);
        %Loop over quantiles of interest
        for q=1:size(Qplot,2)
            patch([(0:H-1) fliplr(0:H-1)],[squeeze(prctile(QIRF_us(k,a,2:end,q,:),alpha/2,[1 2 4 5]))' fliplr(squeeze(prctile(QIRF_us(k,a,2:end,q,:),100-alpha/2,[1 2 4 5]))')],dflt_col{q},'EdgeColor',dflt_col{q},'LineStyle','--','FaceAlpha',.2)
            ymx = max([ymx;max(squeeze(prctile(QIRF_us(k,a,2:end,q,:),100-alpha/2,[1 2 4 5])))]);
            ymn = min([ymn;min(squeeze(prctile(QIRF_us(k,a,2:end,q,:),alpha/2,[1 2 4 5])))]);
            hold on
        end
        yline(0)
        ylm = max(abs(ymx),abs(ymn));
        axis tight; box on;
        hold off
        ylim([-ylm ylm])
        xlim([0 H-1])
        xticks((0:4:H))
        ax = gca();
        Yzr_pos = find(ax.YTick == 0);
        if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
            ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
        end
        if k == K
            xlabel('Quarters ahead','interpreter','latex')
        end
        if k == 1
            if a == gdp_o
                title("Shock to real GDP growth",'interpreter','latex');
            elseif a == irt_o
                title("Shock to Fed. Funds Rate",'interpreter','latex');
            elseif a == inf_o
                title("Shock to CPI inflation",'interpreter','latex');
            elseif a == sri_o
                title("Shock to the fin. cycle",'interpreter','latex');                
            else
                title("Shock to "+ttl{a},'interpreter','latex');
            end
        end
        if a == 1
                if k == irt_o
                    ylabel('Fed. Funds Rate','interpreter','latex');
                elseif k == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif k == inf_o
                    ylabel(["CPI";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{k}],'interpreter','latex');
                end
        end
        set(gca,'TickLabelInterpreter','latex')
    end
end

%Define the legend
lh = legend(...
    round(Quant(Qplot),2) + " quantile (Posterior mean + 95\% credible band)"...
    ,'Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for c = 2:length(fz.Children(1).Children)
    fz.Children.Children(c).YAxis(2).TickValues = [];

    fz.Children.Children(c).FontSize = fz.Children.Children(c).FontSize + 1;
end

fz.Children.Children(1).FontSize = fz.Children.Children(1).FontSize + 1;

print(fig_pth + "\WA_H10",'-dpdf','-vector','-bestfit');

%% END