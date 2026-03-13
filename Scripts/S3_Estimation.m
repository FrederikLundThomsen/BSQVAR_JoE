% ---------------------------------------------------------------- %
% ------------------------- 3. ESTIMATION ------------------------ %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script estimates the US and euro area BSQVARs

% Figures: 
% - E.1
% - G.1-G.9
% - H.1-H.9

% ---------------------------------------------------------------- %

%% ---- Model specification

%Number of lags
P = 4; %Endogenous
S = P; %Exogenous - Note that S=0, i.e. contemporaneous effects from exogenous to endogenous variables, will be present

%Quantile index
Quant = 0.05:0.05:0.95;
nQ = length(Quant);

%Specify the ordering of endogenous variables
ciss_o = 4;
sri_o = 1;
inf_o = 2;
gdp_o = 3;
irt_o = 5;

%Create a variable mapping
var_map = zeros(1,K);
for o = 1:K
    var_map(o) = find([ciss_o sri_o inf_o gdp_o irt_o] == o);
end

%Limits on the time series can be specified here
smpl_lmt_US = US_qe < datetime('2023-01-01','Format','yyyy-MM-dd');
smpl_lmt_EA = qe_dt < datetime('2023-01-01','Format','yyyy-MM-dd');

%Reorder and select the data
Y_US = Yu(smpl_lmt_US,var_map);
Y_EA = Ye(smpl_lmt_EA,var_map);

D_US = D_US(smpl_lmt_US,:);
D_EA = D_EA(smpl_lmt_EA,:);
DM_US = DM_US(smpl_lmt_US,:);
DM_ea = DM_ea(smpl_lmt_EA,:);

X_US = X_US(smpl_lmt_US,:);
X_EA = X_EA(smpl_lmt_EA,:);

qe_dt = qe_dt(smpl_lmt_EA);
US_qe = US_qe(smpl_lmt_US);

%Length of euro area samples
Te = size(Y_EA,1);

%Construct labels for endogenous+exogenous variables
ttla = {'CISS','Fin. cycle','HICP inflation','Real GDP growth','3M OIS rate'};
ttlx = {'Commodity price inflation'};
ttl = ttla(var_map);

%Specify any zero restrictions in the model. These will apply to both the
%US and euro area models as they have the exact same dimensions
ZeroRestrictions = ones(C+K*P+M*(P+1)+K,K);

%Implementation lag of monetary policy on inflation and growth
ZeroRestrictions(C+irt_o,inf_o) = 0;
ZeroRestrictions(C+irt_o,gdp_o) = 0;

%Commodity prices only go through consumer price inflation
ZeroRestrictions(C+K*P+1:M:C+K*P+M*(P+1),~ismember(1:K,inf_o)) = 0;

%% ---- Bayesian modalities

%Number of draws from the posterior distributions
biter = 5e+3;

%Burn-in size
burn = biter/2;

%Thinning value, i.e. the interval between draws from burn to biter used to
%compute relevant statistics, e.g. credible intervals, posterior means etc.
thin = 1;

%Activate looseness parameter
loose = 1;

%% ---- US estimation

%Specify Minnesota prior for all parameters at all quantiles

%Prior means
Dprior = repmat(zeros(K,C).*ZeroRestrictions(1:C,:)',1,1,nQ);
A0prior = zeros(K,K,nQ);

AR_one = [0.9 0.9 1.0 0.9 1.0];
AR_one = AR_one(var_map);

ALprior = repmat([diag(AR_one) zeros(K,K*(P-1))].*ZeroRestrictions(C+1:C+K*P,:)',1,1,nQ);
BLprior = zeros(K,M*(P+1),nQ);

%Pior variance

%General tightness of prior (Small value implies tighter prior)
phi0 = 0.2*1e+0;
%Tightness of prior on other/off-diagonal variables (Small value implies tighter prior)
phi1 = 0.5*1e+0;
%Tightness of prior on exogenous and deterministic variables (Small value implies tighter prior)
phi2 = 1e+5;
%Inverse tightness of lags
phi3 = 1;

%Obtain residuals from standard frequentist quantile autoregression on the median
%to be used for scaling of prior variances

%Initialize vector to store standard deviation of residuals
eps_sd = zeros(1,K);

%Loop over variables
for k = 1:K
    [~,~,~,~,~,~,epsY] = QVAR_system_Z(Y_US(:,k),P,[],0,[D_US(:,1) DM_US],0.5,[]);
    eps_sd(k) = std(epsY);
end

DVprior = repmat(phi0*phi2*ones(K,C).*ZeroRestrictions(1:C,:)',1,1,nQ).^2;
A0Vprior = (phi0*phi1*ones(K,K,nQ).*(eps_sd'*(1./eps_sd))).^2;
ALVprior = repmat(phi0*(phi1*(1-eye(K,K*P)) + eye(K,K*P)).*cell2mat(arrayfun(@(l)(1/(l^phi3))*((eps_sd'*(1./eps_sd))),1:P,'UniformOutput',false)).*ZeroRestrictions(C+1:C+K*P,:)',1,1,nQ).^2;
BLVprior = repmat(phi0*phi2*ones(K,M*(P+1)).*ZeroRestrictions(C+K*P+1:C+K*P+M*(P+1),:)',1,1,nQ).^2;

%Collect the priors in cells
%NOTE: ORDERING IS IMPORTANT! SHOULD BE ORDERED
% i) Deterministic TERMS
% ii) ENDOGENOUS (LAGS WITH COLUMN ORDER [Y1(L1)...YK(L1), Y1(L2)....YK(L2),...YK(LP)])
% iii) EXOGENOUS (SIMULTANEOUS + LAGS WITH COLUMN ORDER [X1(L0)...XM(L0),...XM(LS)])
% iv) ENDOGENOUS (SIMULTANEOUS)

Bprior = {Dprior,ALprior,BLprior,A0prior};
VBprior = {DVprior,ALVprior,BLVprior,A0Vprior};

%Reset RNG
rng(rnd_seed)

%Run the Gibbs sampler, drawing parameters for each variable at each
%quantile, equation-by-equation
[D_us,A0_us,AL_us,BL_us,V_us,eps_us,muUS,sigUS] = QVAR_system_ZB(Y_US,P,X_US,S,D_US,Quant,ZeroRestrictions,biter,Bprior,VBprior,loose);

%Draw from the posterior of the univariate quantile regressions for the
%exogenous variables
b_bbx_US = [];
if ~isempty(X_US)
    DpriorX = repmat(zeros(1,C),1,1,nQ);
    A0priorX = zeros(1,1,nQ);
    ALpriorX = repmat(eye(1,P),1,1,nQ);
    BLpriorX = zeros(1,1,nQ);
    
    DVpriorX = repmat(phi0*phi2*ones(1,C),1,1,nQ).^2;
    A0VpriorX = phi0*phi1*ones(1,1,nQ).^2;
    ALVpriorX = repmat(phi0*(phi1*(1-eye(1,P)) + eye(1,P)).*cell2mat(arrayfun(@(l)(1/(l^phi3)),1:P,'UniformOutput',false)),1,1,nQ).^2;
    BLVpriorX = ones(1,1,nQ).^2;

    BpriorX = {DpriorX,ALpriorX,A0priorX};
    VBpriorX = {DVpriorX,ALVpriorX,A0VpriorX};

    %Loop over exogenous variables
    for m=1:M
        %Reset RNG
        rng(rnd_seed)
        %Draw from the posterior
        [D_bx_us,~,AL_bx_us,~,~,~] = QVAR_system_ZB(X_US(:,m),P,[],0,D_US,Quant,[],biter,BpriorX,VBpriorX,loose);
        %Store the draws in a matrix of size (C+P)*M x biter
        b_bbx_US = [b_bbx_US; reshape(cat(2,D_bx_us,AL_bx_us),[],biter)];
    end
end

%% ---- Euro area estimation

%Use the posterior means and variances from the US estimation as priors for
%euro area estimation

%Take the mean of the relevant US parameters
Dprior = mean(D_us(:,:,:,burn+1:thin:end),length(size(D_us)));
A0prior = mean(A0_us(:,:,:,burn+1:thin:end),length(size(A0_us)));
ALprior = mean(AL_us(:,:,:,burn+1:thin:end),length(size(AL_us)));
if ~isempty(BL_us)
    BLprior = mean(BL_us(:,:,:,burn+1:thin:end),length(size(BL_us)));
else
    BLprior = [];
end

Bprior_ea = {Dprior,ALprior,BLprior,A0prior};
VBprior_ea = {var(D_us(:,:,:,burn+1:thin:end),0,length(size(D_us))),var(AL_us(:,:,:,burn+1:thin:end),0,length(size(AL_us))),var(BL_us(:,:,:,burn+1:thin:end),0,length(size(BL_us))),var(A0_us(:,:,:,burn+1:thin:end),0,length(size(A0_us)))};

%Reset RNG
rng(rnd_seed)

%Posterior draws based on US inference
[D_ea,A0_ea,AL_ea,BL_ea,V_ea,eps_ea,muEA,sigEA] = QVAR_system_ZB(Y_EA,P,X_EA,S,D_EA,Quant,ZeroRestrictions,biter,Bprior_ea,VBprior_ea,loose);

%Draw from the posterior of the univariate quantile regressions for the
%exogenous variables
b_bbx_EA = [];
if ~isempty(X_EA)
    DpriorX = mean(D_bx_us(:,:,:,burn+1:thin:end),length(size(D_bx_us)));
    A0priorX = zeros(1,1,nQ);
    ALpriorX = mean(AL_bx_us(:,:,:,burn+1:thin:end),length(size(D_bx_us)));
    BLpriorX = zeros(1,1,nQ);
    
    DVpriorX = var(D_bx_us(:,:,:,burn+1:thin:end),0,length(size(D_bx_us)));
    A0VpriorX = phi0*phi1*ones(1,1,nQ).^2;
    ALVpriorX = var(AL_bx_us(:,:,:,burn+1:thin:end),0,length(size(D_bx_us)));
    BLVpriorX = ones(1,1,nQ).^2;

    BpriorX = {DpriorX,ALpriorX,A0priorX};
    VBpriorX = {DVpriorX,ALVpriorX,A0VpriorX};

    for m=1:M
        %Reset RNG
        rng(rnd_seed)
        %Draw from the posterior
        [D_bx_ea,~,AL_bx_ea,~,~,~] = QVAR_system_ZB(X_EA(:,m),P,[],0,D_EA,Quant,[],biter,BpriorX,VBpriorX,loose);
        %Store the draws in a matrix of size (C+P)*M x biter
        b_bbx_EA = [b_bbx_EA; reshape(cat(2,D_bx_ea,AL_bx_ea),[],biter)];
    end
end

%% ---- Figures
%% ---- ---- Web appendix

%Credible level to plot in percentages
alpha = 5;

%Estimate the linear VAR version of the US model
[~,D_us_l,A0_us_l,AL_us_l,BL_us_l,~,~] = SVAR_system_Z(Y_US,P,X_US,S,D_US,ZeroRestrictions);

%Estimate the linear VAR version of the euro area model
[~,D_ea_l,A0_ea_l,AL_ea_l,BL_ea_l,~,~] = SVAR_system_Z(Y_EA,P,X_EA,S,D_EA,ZeroRestrictions);

%% ---- ---- Figure E.1: Euro area and U.S. time series

%Call a new figure
figure
tiledlayout(6,2,"TileSpacing","compact"',"Padding","compact")

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1;
figp(4) = 21*3/5*2;
set(gcf,'Position',figp);    

%Plot the individual time series
nexttile(1)
plot(qe_dt(smpl_lmt_EA),Y_EA(:,sri_o),'LineWidth',2)
axis tight
box on
title('Euro area','interpreter','latex')
ylabel('Financial cycle','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(2)
plot(US_qe(smpl_lmt_US),Y_US(:,sri_o),'LineWidth',2)
axis tight
box on
title('U.S.','interpreter','latex')
ylabel('Financial cycle','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(3)
plot(qe_dt(smpl_lmt_EA),Y_EA(:,inf_o),'LineWidth',2)
axis tight
box on
ylabel({'HICP inflation' ; '\% annualised'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(4)
plot(US_qe(smpl_lmt_US),Y_US(:,inf_o),'LineWidth',2)
axis tight
box on
ylabel({'CPI inflation' ; '\% annualised'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(5)
plot(qe_dt(smpl_lmt_EA),Y_EA(:,gdp_o),'LineWidth',2)
axis tight
box on
ylim([-10 10])
ylabel({'Real GDP growth' ; '\% annualised'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(6)
plot(US_qe(smpl_lmt_US),Y_US(:,gdp_o),'LineWidth',2)
axis tight
box on
ylim([-10 10])
ylabel({'Real GDP growth' ; '\% annualised'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(7)
plot(qe_dt(smpl_lmt_EA),Y_EA(:,ciss_o),'LineWidth',2)
axis tight
box on
ylim([0 1])
ylabel({'CISS'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(8)
plot(US_qe(smpl_lmt_US),Y_US(:,ciss_o),'LineWidth',2)
axis tight
box on
ylim([0 1])
ylabel({'CISS'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')    
yline(0)

nexttile(9)
plot(qe_dt(smpl_lmt_EA),Y_EA(:,irt_o),'LineWidth',2)
axis tight
box on
ylabel({'3-month OIS rate' ; '\% p.a.'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
yline(0)

nexttile(10)
plot(US_qe(smpl_lmt_US),Y_US(:,irt_o),'LineWidth',2)
axis tight
box on
ylabel({'Federal Funds Rate' ; '\% p.a.'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')        
yline(0)

nexttile(11,[1 2])
plot(US_qe(smpl_lmt_US),X_US(:,1),'LineWidth',2)
axis tight
box on
ylabel({'Cmdt. price inflation' ; '\% annualised'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')        
yline(0)

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 1:length(fz.Children(1).Children)
    fz.Children(1).Children(f).FontSize = fz.Children(1).Children(f).FontSize - 1;
    fz.Children(1).Children(f).XAxis(1).Label.FontSize = fz.Children(1).Children(f).XAxis(1).Label.FontSize +1;
    fz.Children(1).Children(f).Title.FontSize = fz.Children(1).Children(f).Title.FontSize +1;
    fz.Children(1).Children(f).YAxis(2).TickValues = [];

end

%Print
print(fig_pth+"\WA_E1",'-dpdf','-vector','-bestfit') ;

%% ---- ---- Figure G.1: Posterior inference for the euro area, omega(gamma)

%Call a new figure
figure
tiledlayout(1,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.5;
figp(4) = 21*3/5*1;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','portrait')

%Loop over endogenous variables
for a=1:K
    for b = 1:1
        nexttile
        plot(Quant,mean(squeeze(D_ea(a,b,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
        hold on
        patch([Quant fliplr(Quant)],[prctile(squeeze(D_ea(a,b,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(D_ea(a,b,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
        hold on
        plot(Quant,ones(nQ,1)*D_ea_l(a,b),'-.','LineWidth',2,'Color',dflt_col{2})
        yline(0)
        axis tight; box on;
        hold off
        xlim([Quant(1) Quant(end)])
        ax = gca();
        Yzr_pos = find(ax.YTick == 0);
        if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
            ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
        end
        xlabel('Quantile','interpreter','latex');
        if b == 1
            if a == irt_o
                title('3M OIS rate','interpreter','latex');
            elseif a == gdp_o
                title(["Real GDP";"growth"],'interpreter','latex');
            elseif a == inf_o
                title(["HICP";"inflation"],'interpreter','latex');
            else
                title([ttl{a}],'interpreter','latex');
            end
        end
    end
end

set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend('Posterior mean',100-alpha + "\% credible interval",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 2:length(fz.Children(1).Children)
    fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
    fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize + 1;
    fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize + 1;
    fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize + 1;
    
    fz.Children.Children(f).YAxis(2).TickValues = [];
end

for c = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(c).Type,'axes')
        yyaxis(fz.Children(1).Children(c),'left')
        if ~isempty(fz.Children(1).Children(c).Children)
            fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
        end
    end
end

%Print
print(fig_pth+"\WA_G1",'-dpdf','-vector','-bestfit');

%% ---- ---- Figure G.2: Posterior inference for the euro area, A0(gamma)

%Call a new figure
figure
tiledlayout(K,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.5;
figp(4) = 21*3/5*1.5;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','portrait')

%Loop over endogenous variables
for a=1:K
    %Loop over endogenous variables on the RHS
    for b=1:K
        if a > b
            nexttile((a-1)*K + b)
            plot(Quant,mean(squeeze(A0_ea(a,b,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
            hold on
            patch([Quant fliplr(Quant)],[prctile(squeeze(A0_ea(a,b,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(A0_ea(a,b,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
            hold on
            plot(Quant,ones(nQ,1)*A0_ea_l(a,b),'-.','LineWidth',2,'Color',dflt_col{2})
            yline(0);
            axis tight; box on;
            hold off
            xlim([Quant(1) Quant(end)])
            ax = gca();
            Yzr_pos = find(ax.YTick == 0);
            if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
                ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
            end
            if a == b+1
                title([ttl{b}],'interpreter','latex');
            end
            if a == K
                xlabel('Quantile','interpreter','latex');
            end
            if b == 1
                if a == irt_o
                    ylabel('3M OIS rate','interpreter','latex');
                elseif a == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif a == inf_o
                    ylabel(["HICP";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{a}],'interpreter','latex');
                end
            end
        end
    end
end

set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend('Posterior mean',100-alpha + "\% credible interval",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 2:length(fz.Children(1).Children)
    fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
    fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
    fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
    fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
    fz.Children.Children(f).YAxis(2).TickValues = [];
end

for c = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(c).Type,'axes')
        yyaxis(fz.Children(1).Children(c),'left')
        if ~isempty(fz.Children(1).Children(c).Children)
            fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
        end
    end
end

%Print
print(fig_pth+"\WA_G2",'-dpdf','-vector','-bestfit');

%% ---- ---- Figure G.3-G.6: Posterior inference for the euro area, A1(gamma)-A4(gamma)

%Loop over lags
for p=1:P
    %Call new figure
    figure
    tiledlayout(K,K,'TileSpacing','Compact','Padding','Compact');

    %Scaling
    set(gcf,'Units','centimeters')
    figp = get(gcf,'Position');
    figp(3) = 21*1.5;
    figp(4) = 21*3/5*1.5;
    set(gcf,'Position',figp);
    set(gcf,'PaperOrientation','portrait')

    %Loop over endogenous variables
    for a=1:K
        %Loop over endogenous variables on the RHS
        for b=1:K
            nexttile
            plot(Quant,mean(squeeze(AL_ea(a,b+(p-1)*K,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
            hold on
            patch([Quant fliplr(Quant)],[prctile(squeeze(AL_ea(a,b+(p-1)*K,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(AL_ea(a,b+(p-1)*K,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
            hold on
            plot(Quant,ones(nQ,1)*AL_ea_l(a,b+(p-1)*K),'-.','LineWidth',2,'Color',dflt_col{2})
            yline(0)
            axis tight; box on;
            hold off
            xlim([Quant(1) Quant(end)])
            ax = gca();
            Yzr_pos = find(ax.YTick == 0);
            if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
                ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
            end
            if a == 1
                if b == K
                    title('3M OIS rate','interpreter','latex');
                else
                    title([ttl{b}],'interpreter','latex');
                end
            end
            if a == K
                xlabel('Quantile','interpreter','latex');
            end
            if b == 1
                if a == irt_o
                    ylabel('3M OIS rate','interpreter','latex');
                elseif a == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif a == inf_o
                    ylabel(["HICP";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{a}],'interpreter','latex');
                end
            end
        end
    end
    set(gca,'TickLabelInterpreter','latex')

    %Define the legend
    lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
    lh.Layout.Tile = 'North';

    %Call function to format active figure
    PlotFmt([],1)

    %Ad hoc modifications
    fz = gcf();
    for f = 2:length(fz.Children(1).Children)
        fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
        fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
        fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
        fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
        fz.Children.Children(f).YAxis(2).TickValues = [];
    end

    for c = 1:length(fz.Children(1).Children)
        if strcmp(fz.Children(1).Children(c).Type,'axes')
            yyaxis(fz.Children(1).Children(c),'left')
            if ~isempty(fz.Children(1).Children(c).Children)
                fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
            end
        end
    end
    
    %Print
    print(fig_pth+"\WA_G"+(p+2),'-dpdf','-vector','-bestfit');

end

%% ---- ---- Figure G.7: Posterior inference for the euro area exogenous variables, C(gamma)

if ~isempty(X_EA)

    %Call a new figure
    figure
    tiledlayout(K,S+1,'TileSpacing','Compact','Padding','Compact');

    %Scaling
    set(gcf,'Units','centimeters')
    figp = get(gcf,'Position');
    figp(3) = 21*1.5;
    figp(4) = 21*3/5*1.5;
    set(gcf,'Position',figp);
    set(gcf,'PaperOrientation','portrait')

    %Loop over endogenous variables
    for a=1:K
        %Loop over lags of exogenous variables
        for b=0:S
            nexttile
            plot(Quant,mean(squeeze(BL_ea(a,b+1,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
            hold on
            patch([Quant fliplr(Quant)],[prctile(squeeze(BL_ea(a,b+1,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(BL_ea(a,b+1,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
            hold on
            plot(Quant,ones(nQ,1)*BL_ea_l(a,b+1),'-.','LineWidth',2,'Color',dflt_col{2})
            yline(0)
            axis tight; box on;
            hold off
            xlim([Quant(1) Quant(end)])
            ax = gca();
            Yzr_pos = find(ax.YTick == 0);
            if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
                ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
            end
            if a == 1
                title("Commodity prices: Lag " + b,'interpreter','latex');
            end
            if a == K
                xlabel('Quantile','interpreter','latex');
            end
            if b == 0
                if a == irt_o
                    ylabel('3M OIS rate','interpreter','latex');
                elseif a == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif a == inf_o
                    ylabel(["HICP";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{a}],'interpreter','latex');
                end
            end
        end
    end
    set(gca,'TickLabelInterpreter','latex')

    %Define the legend
    lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
    lh.Layout.Tile = 'North';

    %Call function to format active figure
    PlotFmt([],1)

    %Ad hoc modifications
    fz = gcf();
    for f = 2:length(fz.Children(1).Children)
        fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
        fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
        fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
        fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
       
        fz.Children.Children(f).YAxis(2).TickValues = [];
    end

    for c = 1:length(fz.Children(1).Children)
        if strcmp(fz.Children(1).Children(c).Type,'axes')
            yyaxis(fz.Children(1).Children(c),'left')
            if ~isempty(fz.Children(1).Children(c).Children)
                fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
            end
        end
    end

    %Print
    print(fig_pth+"\WA_G7",'-dpdf','-vector','-bestfit');

end

%% ---- ---- Figure G.8: Posterior inference for the euro area dummy variables, B(gamma)

%Call a new figure
figure
tiledlayout(C-1,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.5;
figp(4) = 21*3/5*1.5;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','portrait')

%Loop over deterministic terms (sans the constant)
for b=2:C
    %Loop over endogeous variables
    for a = 1:K
        nexttile
        plot(Quant,mean(squeeze(D_ea(a,b,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
        hold on
        patch([Quant fliplr(Quant)],[prctile(squeeze(D_ea(a,b,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(D_ea(a,b,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
        hold on
        plot(Quant,ones(nQ,1)*D_ea_l(a,b),'-.','LineWidth',2,'Color',dflt_col{2})
        yline(0)
        axis tight; box on;
        hold off
        xlim([Quant(1) Quant(end)])
        ax = gca();
        Yzr_pos = find(ax.YTick == 0);
        if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
            ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
        end
        if b == C
            xlabel('Quantile','interpreter','latex');
        end
        if a == 1
            if b == 1
                ylabel('Constant','interpreter','latex');
            elseif b == 2
                ylabel('2020Q1 dummy','interpreter','latex');
            elseif b == 3
                ylabel('2020Q2 dummy','interpreter','latex');
            elseif b == 4
                ylabel('2020Q3 dummy','interpreter','latex');
            elseif b == 5
                ylabel('2020Q4 dummy','interpreter','latex');
            end
        end           
        if b == 2
            if a == irt_o
                title('3M OIS rate','interpreter','latex');
            elseif a == gdp_o
                title(["Real GDP";"growth"],'interpreter','latex');
            elseif a == inf_o
                title(["HICP";"inflation"],'interpreter','latex');
            else
                title([ttl{a}],'interpreter','latex');
            end
        end
    end
end

set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 2:length(fz.Children(1).Children)
    fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
    fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
    fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
    fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
    fz.Children.Children(f).YAxis(2).TickValues = [];
end

for c = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(c).Type,'axes')
        yyaxis(fz.Children(1).Children(c),'left')
        if ~isempty(fz.Children(1).Children(c).Children)
            fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
        end
    end
end

%Print
print(fig_pth+"\WA_G8",'-dpdf','-vector','-bestfit');

%% ---- ---- Figure G.9: Posterior means for the euro area tightness parameter, lambda(gamma)

%Call a new figure
figure
tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 25;
figp(4) = 15;
set(gcf,'Position',figp);

%Plot the posterior means of the looseness paramamter, 
nexttile
hold on
plot(Quant,mean(squeeze(muEA(sri_o,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
plot(Quant,mean(squeeze(muEA(inf_o,:,burn+1:thin:end)),2),'--','LineWidth',2,'Color',dflt_col{2})
plot(Quant,mean(squeeze(muEA(gdp_o,:,burn+1:thin:end)),2),':','LineWidth',2,'Color',dflt_col{3})
plot(Quant,mean(squeeze(muEA(ciss_o,:,burn+1:thin:end)),2),'-o','LineWidth',2,'Color',dflt_col{4},'MarkerFaceColor',dflt_col{4})
plot(Quant,mean(squeeze(muEA(irt_o,:,burn+1:thin:end)),2),'-d','LineWidth',2,'Color',dflt_col{5},'MarkerFaceColor',dflt_col{5})
hold off
yline(0);
axis tight; box on;
hold off
xlim([Quant(1) Quant(end)])
xlabel('Quantile','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')        

%Define the legend
lh = legend('SRI','HICP inflation','Real GDP growth','CISS','3-month OIS rate','Orientation','Horizontal','Location','Northoutside','NumColumns',K,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();

fz.Children(1).Children(2).YAxis(2).TickValues = [];

%Print
print(fig_pth+"\WA_G9",'-dpdf','-vector','-bestfit') ;

%% ---- ---- Figure H.1: Posterior inference for the U.S., omega(gamma)

%Call a new figure
figure
tiledlayout(1,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.5;
figp(4) = 21*3/5*1;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','portrait')

%Loop over endogenous variables
for a=1:K
    for b = 1:1
        nexttile
        plot(Quant,mean(squeeze(D_us(a,b,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
        hold on
        patch([Quant fliplr(Quant)],[prctile(squeeze(D_us(a,b,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(D_us(a,b,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
        hold on
        plot(Quant,ones(nQ,1)*D_us_l(a,b),'-.','LineWidth',2,'Color',dflt_col{2})
        yline(0)
        axis tight; box on;
        hold off
        xlim([Quant(1) Quant(end)])
        ax = gca();
        Yzr_pos = find(ax.YTick == 0);
        if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
            ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
        end
        xlabel('Quantile','interpreter','latex');
        if b == 1
            if a == irt_o
                title('Fed. Funds Rate','interpreter','latex');
            elseif a == gdp_o
                title(["Real GDP";"growth"],'interpreter','latex');
            elseif a == inf_o
                title(["CPI";"inflation"],'interpreter','latex');
            else
                title([ttl{a}],'interpreter','latex');
            end
        end
    end
end

set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 2:length(fz.Children(1).Children)
    fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
    fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
    fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
    fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
    fz.Children.Children(f).YAxis(2).TickValues = [];
end

for c = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(c).Type,'axes')
        yyaxis(fz.Children(1).Children(c),'left')
        if ~isempty(fz.Children(1).Children(c).Children)
            fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
        end
    end
end

%Print
print(fig_pth+"\WA_H1",'-dpdf','-vector','-bestfit');

%% ---- ---- Figure H.2: Posterior inference for the U.S., A0(gamma)

%Call new figure
figure
tiledlayout(K,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.5;
figp(4) = 21*3/5*1.5;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','portrait')

%Loop over endogenous variables
for a=1:K
    %Loop over endogenous variables on the RHS
    for b=1:K
        if a > b
            nexttile((a-1)*K + b)
            plot(Quant,mean(squeeze(A0_us(a,b,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
            hold on
            patch([Quant fliplr(Quant)],[prctile(squeeze(A0_us(a,b,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(A0_us(a,b,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
            hold on
            plot(Quant,ones(nQ,1)*A0_us_l(a,b),'-.','LineWidth',2,'Color',dflt_col{2})
            yline(0)
            axis tight; box on;
            hold off
            xlim([Quant(1) Quant(end)])
            ax = gca();
            Yzr_pos = find(ax.YTick == 0);
            if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
                ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
            end
            if a == b+1
                title([ttl{b}],'interpreter','latex');
            end
            if a == K
                xlabel('Quantile','interpreter','latex');
            end
            if b == 1
                if a == irt_o
                    ylabel('Fed. Funds Rate','interpreter','latex');
                elseif a == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif a == inf_o
                    ylabel(["CPI";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{a}],'interpreter','latex');
                end
            end
        end
    end
end

set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 2:length(fz.Children(1).Children)
    fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
    fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
    fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
    fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
    fz.Children.Children(f).YAxis(2).TickValues = [];
end

for c = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(c).Type,'axes')
        yyaxis(fz.Children(1).Children(c),'left')
        if ~isempty(fz.Children(1).Children(c).Children)
            fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
        end
    end
end

%Print
print(fig_pth+"\WA_H2",'-dpdf','-vector','-bestfit');

%% ---- ---- Figure H.3-H.6: Posterior inference for the U.S., A1(gamma)-A4(gamma)

%Loop over lags
for p=1:P
    %Call a new figure
    figure
    tiledlayout(K,K,'TileSpacing','Compact','Padding','Compact');

    %Scaling
    set(gcf,'Units','centimeters')
    figp = get(gcf,'Position');
    figp(3) = 21*1.5;
    figp(4) = 21*3/5*1.5;
    set(gcf,'Position',figp);
    set(gcf,'PaperOrientation','portrait')

    %Loop over endogenous variables
    for a=1:K
        %Loop over endogenous variables on the RHS
        for b=1:K
            nexttile
            plot(Quant,mean(squeeze(AL_us(a,b+(p-1)*K,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
            hold on
            patch([Quant fliplr(Quant)],[prctile(squeeze(AL_us(a,b+(p-1)*K,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(AL_us(a,b+(p-1)*K,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
            hold on
            plot(Quant,ones(nQ,1)*AL_us_l(a,b+(p-1)*K),'-.','LineWidth',2,'Color',dflt_col{2})
            yline(0)
            axis tight; box on;
            hold off
            xlim([Quant(1) Quant(end)])
            ax = gca();
            Yzr_pos = find(ax.YTick == 0);
            if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
                ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
            end
            if a == 1
                if b == irt_o
                    title('Fed. Funds Rate','interpreter','latex');
                elseif b == inf_o
                    title('CPI inflation','interpreter','latex');
                else
                    title([ttl{b}],'interpreter','latex');
                end
            end
            if a == K
                xlabel('Quantile','interpreter','latex');
            end
            if b == 1
                if a == irt_o
                    ylabel('Fed. Funds Rate','interpreter','latex');
                elseif a == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif a == inf_o
                    ylabel(["CPI";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{a}],'interpreter','latex');
                end
            end
        end
    end
    set(gca,'TickLabelInterpreter','latex')

    %Define the legend
    lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
    lh.Layout.Tile = 'North';

    %Call function to format active figure
    PlotFmt([],1)

    %Ad hoc modifications
    fz = gcf();
    for f = 2:length(fz.Children(1).Children)
        fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
        fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
        fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
        fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
        fz.Children.Children(f).YAxis(2).TickValues = [];
    end

    for c = 1:length(fz.Children(1).Children)
        if strcmp(fz.Children(1).Children(c).Type,'axes')
            yyaxis(fz.Children(1).Children(c),'left')
            if ~isempty(fz.Children(1).Children(c).Children)
                fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
            end
        end
    end
        
    %Print
    print(fig_pth+"\WA_H"+(p+2),'-dpdf','-vector','-bestfit');

end

%% ---- ---- Figure H.7: Posterior inference for the U.S. exogenous variables, C(gamma)

if ~isempty(X_US)

    %Call a new figure
    figure
    tiledlayout(K,S+1,'TileSpacing','Compact','Padding','Compact');

    %Scaling
    set(gcf,'Units','centimeters')
    figp = get(gcf,'Position');
    figp(3) = 21*1.5;
    figp(4) = 21*3/5*1.5;
    set(gcf,'Position',figp);
    set(gcf,'PaperOrientation','portrait')

    %Loop over endogenous variables
    for a=1:K
        %Loop over lags of exogenous variables
        for b=0:S
            nexttile
            plot(Quant,mean(squeeze(BL_us(a,b+1,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
            hold on
            patch([Quant fliplr(Quant)],[prctile(squeeze(BL_us(a,b+1,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(BL_us(a,b+1,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
            hold on
            plot(Quant,ones(nQ,1)*BL_us_l(a,b+1),'-.','LineWidth',2,'Color',dflt_col{2})
            yline(0)
            axis tight; box on;
            hold off
            xlim([Quant(1) Quant(end)])
            ax = gca();
            Yzr_pos = find(ax.YTick == 0);
            if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
                ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
            end
            if a == 1
                title("Commodity prices: Lag " + b,'interpreter','latex');
            end
            if a == K
                xlabel('Quantile','interpreter','latex');
            end
            if b == 0
                if a == irt_o
                    ylabel('Fed. Funds Rate','interpreter','latex');
                elseif a == gdp_o
                    ylabel(["Real GDP";"growth"],'interpreter','latex');
                elseif a == inf_o
                    ylabel(["CPI";"inflation"],'interpreter','latex');
                else
                    ylabel([ttl{a}],'interpreter','latex');
                end
            end
        end
    end
    set(gca,'TickLabelInterpreter','latex')

    %Define the legend
    lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
    lh.Layout.Tile = 'North';

    %Call function to format active figure
    PlotFmt([],1)

    %Ad hoc modifications
    fz = gcf();
    for f = 2:length(fz.Children(1).Children)
        fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
        fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
        fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
        fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
        fz.Children.Children(f).YAxis(2).TickValues = [];
    end

    for c = 1:length(fz.Children(1).Children)
        if strcmp(fz.Children(1).Children(c).Type,'axes')
            yyaxis(fz.Children(1).Children(c),'left')
            if ~isempty(fz.Children(1).Children(c).Children)
                fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
            end
        end
    end
    
    %Print
    print(fig_pth+"\WA_H7",'-dpdf','-vector','-bestfit');
end

%% ---- ---- Figure H.8: Posterior inference for the U.S. dummy variables, B(gamma)

%Call a new figure
figure
tiledlayout(C-1,K,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.5;
figp(4) = 21*3/5*1.5;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','portrait')

%Loop over deterministic terms (sans the constant)
for b=2:C
    %Loop over endogenous variables
    for a = 1:K
        nexttile
        plot(Quant,mean(squeeze(D_us(a,b,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
        hold on
        patch([Quant fliplr(Quant)],[prctile(squeeze(D_us(a,b,:,burn+1:thin:end)),alpha/2,2)' flip(prctile(squeeze(D_us(a,b,:,burn+1:thin:end)),100-alpha/2,2))'],'--','LineWidth',2,'FaceColor',dflt_col{1},'EdgeColor',dflt_col{1},'LineStyle','--','FaceAlpha',.2,'LineWidth',1);
        hold on
        plot(Quant,ones(nQ,1)*D_us_l(a,b),'-.','LineWidth',2,'Color',dflt_col{2})
        yline(0)
        axis tight; box on;
        hold off
        xlim([Quant(1) Quant(end)])
        ax = gca();
        Yzr_pos = find(ax.YTick == 0);
        if (length(ax.YTick) > 4)||((length(ax.YTick) > 3)&&((min(ax.YTick) == 0)||(max(ax.YTick == 0))))
            ax.YTick = ax.YTick(2-mod(Yzr_pos,2):2:end);                
        end
        if b == C
            xlabel('Quantile','interpreter','latex');
        end
        if a == 1
            if b == 1
                ylabel('Constant','interpreter','latex');
            elseif b == 2
                ylabel('2020Q1 dummy','interpreter','latex');
            elseif b == 3
                ylabel('2020Q2 dummy','interpreter','latex');
            elseif b == 4
                ylabel('2020Q3 dummy','interpreter','latex');
            elseif b == 5
                ylabel('2020Q4 dummy','interpreter','latex');
            end
        end           
        if b == 2
            if a == irt_o
                title('Fed. Funds Rate','interpreter','latex');
            elseif a == gdp_o
                title(["Real GDP";"growth"],'interpreter','latex');
            elseif a == inf_o
                title(["CPI";"inflation"],'interpreter','latex');
            else
                title([ttl{a}],'interpreter','latex');
            end
        end
    end
end

set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend('Posterior mean',100-alpha + "\% credible band",'Least squares estimate','Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
for f = 2:length(fz.Children(1).Children)
    fz.Children.Children(f).FontSize = fz.Children.Children(f).FontSize - 2;
    fz.Children.Children(f).YAxis(1).Label.FontSize = fz.Children.Children(f).YAxis(1).Label.FontSize +1;
    fz.Children.Children(f).XAxis(1).Label.FontSize = fz.Children.Children(f).XAxis(1).Label.FontSize +1;
    fz.Children.Children(f).Title.FontSize = fz.Children.Children(f).Title.FontSize +1;
    
    fz.Children.Children(f).YAxis(2).TickValues = [];
end

for c = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(c).Type,'axes')
        yyaxis(fz.Children(1).Children(c),'left')
        if ~isempty(fz.Children(1).Children(c).Children)
            fz.Children(1).Children(c).Children(2).Color = dflt_col{3};
        end
    end
end

%Print
print(fig_pth+"\WA_H8",'-dpdf','-vector','-bestfit');

%% ---- ---- Figure H.9: Posterior mean estimates for U.S. prior tightness parameter i(gamma)

%Call a new figure
figure
tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 25;
figp(4) = 15;
set(gcf,'Position',figp);

%Plot the posterior means of the looseness paramamter, 
nexttile
hold on
plot(Quant,mean(squeeze(muUS(sri_o,:,burn+1:thin:end)),2),'-','LineWidth',2,'Color',dflt_col{1})
plot(Quant,mean(squeeze(muUS(inf_o,:,burn+1:thin:end)),2),'--','LineWidth',2,'Color',dflt_col{2})
plot(Quant,mean(squeeze(muUS(gdp_o,:,burn+1:thin:end)),2),':','LineWidth',2,'Color',dflt_col{3})
plot(Quant,mean(squeeze(muUS(ciss_o,:,burn+1:thin:end)),2),'-o','LineWidth',2,'Color',dflt_col{4},'MarkerFaceColor',dflt_col{4})
plot(Quant,mean(squeeze(muUS(irt_o,:,burn+1:thin:end)),2),'-d','LineWidth',2,'Color',dflt_col{5},'MarkerFaceColor',dflt_col{5})
hold off
yline(0);
axis tight; box on;
hold off
xlim([Quant(1) Quant(end)])
xlabel('Quantile','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')        

%Define the legend
lh = legend('SRI','CPI inflation','Real GDP growth','CISS','Federal Funds Rate','Orientation','Horizontal','Location','Northoutside','NumColumns',K,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();

fz.Children(1).Children(2).YAxis(2).TickValues = [];

%Print
print(fig_pth+"\WA_H9",'-dpdf','-vector','-bestfit') ;    
%% END