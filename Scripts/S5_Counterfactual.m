% ---------------------------------------------------------------- %
% ------------------- 5. Counterfactual policy ------------------- %
% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

% This script produces posterior inference for the counterfactual outcomes
% under different macroprudential policy settings

% Figures: 
% - 2
% - G.10

% ---------------------------------------------------------------- %

%% ---- Counterfactual macroprudential policy

%Number of posterior draws
pDraw = 10e+3;

%Array to store realised quantiles
quant_real = nan(Te-P,K,pDraw);
%Counterfactual values for 'aggressive' policy
y_count_med = nan(Te,K,pDraw);
%Counterfactual values for the 'moderate' policy
y_count_s_int = y_count_med;

%Pick the policy adjustment to close the gap between the median and 
% realized quantiles
kappa_q = 0.5;

%Start a timer
TIME = tic;

%Draw
for bp = 1:pDraw

    %Print progress
    "[" + string(datetime()) + "] " + "Time elapsed: " + sprintf("%s",duration([0, 0, toc(TIME)])) + " Iteration: " + bp + "."
    
    %Set random seed
    rng(rnd_seed + bp)
    pstrr = randi([burn+1 biter],1,1);
    
    %Draw posterior coefficients
    D_bb = D_ea(:,:,:,pstrr);
    A0_bb = A0_ea(:,:,:,pstrr);
    AL_bb = AL_ea(:,:,:,pstrr);
    if M > 0
        BL_bb = BL_ea(:,:,:,pstrr);
        b_bb = b_bbx_EA(:,pstrr);
    else
        BL_bb = [];
        b_bb = [];
    end
    
    %Initialize quantile realizations conditional on full sample estimates
    quant_real_p = nan(Te-P,K);
    %Loop over time periods
    for t = P+1:Te
        %Loop over variables
        for k = 1:K
            %Array to store quantile predictions
            ykp = nan(nQ,1);
            %Loop over quantiles
            for q = 1:nQ
                %Compute quantile prediction one step ahead
                if M > 0
                    ykp(q) = squeeze(D_bb(k,:,q))*D_EA(t,:)' + [squeeze(A0_bb(k,:,q)) squeeze(AL_bb(k,:,q))]*reshape(flip(Y_EA(t-P:t,:))',[],1) + squeeze(BL_bb(k,:,q))*reshape(flip(X_EA(t-P:t,:))',[],1);
                else
                    ykp(q) = squeeze(D_bb(k,:,q))*D_EA(t,:)' + [squeeze(A0_bb(k,:,q)) squeeze(AL_bb(k,:,q))]*reshape(flip(Y_EA(t-P:t,:))',[],1);
                end
            end
            %Compute residual
            ekp = Y_EA(t,k) - ykp;
            %Find the smallest absolute value residual
            ekp_mn = abs(ekp) == min(abs(ekp));
            ekp_mn_num = find(ekp_mn);
            %If than one quantile provides the smallest residual, pick one
            %randomly
            if numel(ekp_mn_num) > 1
                ekp_rnd = randsample(numel(ekp_mn_num),1);
                ekp_pck = ekp_mn_num(ekp_rnd);
            else
                ekp_pck = ekp_mn_num;
            end
            %Store the pick
            quant_real_p(t-P,k) = ekp_pck;
        end
    end
        
    %Now the counterfactual for median realizations
    quant_med_count = quant_real_p;
    quant_med_count(:,sri_o) = find(round(Quant,2) == 0.50); 
    y_count_med_p = [Y_EA(1:P,:) ; zeros(Te-P,K)];
    for t = P+1:Te
       %Loop over variables
        for k = 1:K
            qnt_t = quant_med_count(t-P,k);
            if M > 0
                y_count_med_p(t,k) = squeeze(D_bb(k,:,qnt_t))*D_EA(t,:)' + [squeeze(A0_bb(k,:,qnt_t)) squeeze(AL_bb(k,:,qnt_t))]*reshape(flip(y_count_med_p(t-P:t,:))',[],1) + squeeze(BL_bb(k,:,qnt_t))*reshape(flip(X_EA(t-P:t,:))',[],1);
            else
                y_count_med_p(t,k) = squeeze(D_bb(k,:,qnt_t))*D_EA(t,:)' + [squeeze(A0_bb(k,:,qnt_t)) squeeze(AL_bb(k,:,qnt_t))]*reshape(flip(y_count_med_p(t-P:t,:))',[],1);
            end
        end
    end

    %Now the counterfactual for weighted quantile realizations
    quant_s_int_count = quant_real_p;
    %Pick a weighted quantile
    quant_s_int_count(:,sri_o) = kappa_q*quant_real_p(:,sri_o) + (1-kappa_q)*quant_med_count(:,sri_o); 
    %Randomly round up or down when need be
    quant_s_int_count(:,sri_o) = floor(quant_s_int_count(:,sri_o) + (mod(quant_s_int_count(:,sri_o),1) ~= 0).*(randsample([0 1],Te-P,true)' - mod(quant_s_int_count(:,sri_o),1)));
    %Take care of the edges
    quant_s_int_count(:,sri_o) = max(1,min(nQ,quant_s_int_count(:,sri_o)));
    y_count_s_int_p = y_count_med_p;
    for t = P+1:Te
       %Loop over variables
        for k = 1:K
            qnt_t = quant_s_int_count(t-P,k);
            if M > 0
                y_count_s_int_p(t,k) = squeeze(D_bb(k,:,qnt_t))*D_EA(t,:)' + [squeeze(A0_bb(k,:,qnt_t)) squeeze(AL_bb(k,:,qnt_t))]*reshape(flip(y_count_s_int_p(t-P:t,:))',[],1) + squeeze(BL_bb(k,:,qnt_t))*reshape(flip(X_EA(t-P:t,:))',[],1);
            else
                y_count_s_int_p(t,k) = squeeze(D_bb(k,:,qnt_t))*D_EA(t,:)' + [squeeze(A0_bb(k,:,qnt_t)) squeeze(AL_bb(k,:,qnt_t))]*reshape(flip(y_count_s_int_p(t-P:t,:))',[],1);
            end
        end
    end

    %Store
    quant_real(:,:,bp) = quant_real_p;
    y_count_med(:,:,bp) = y_count_med_p;
    y_count_s_int(:,:,bp) = y_count_s_int_p;

end

%% ---- Figures

%Credible level
alphau = 32;

%% ---- Figure 2: The financial cycle in two different experiments

%Scenario labels
scen_lbl = {'``Aggressive" policy','``Moderate" policy'};

var_plt = sri_o;

%Call a new figure
figure
tiledlayout(1,1,"TileSpacing","compact"',"Padding","compact")

%Scaling
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 21*1.25;
figp(4) = 21*3/5*1*1;
set(gcf,'Position',figp);

%Begin plotting
nexttile
hold on

%Plot realized values
plot(qe_dt,Y_EA(:,var_plt),'-','Color',dflt_col{1},'LineWidth',1.5)

%Plot counterfactual posterior means
for s = 1:2
    y_scen = y_count_med*(s == 1) + y_count_s_int*(s == 2);
    plot(qe_dt,mean(y_scen(:,var_plt,:),[2 3]),'-','Color',dflt_col{s+1},'LineWidth',1.5)
end

%Plot credible bands for the 'aggressive' policy
y_scen = y_count_med;
patch([qe_dt ; flip(qe_dt)],[prctile(y_scen(:,var_plt,:),alphau/2,[2 3]) ; flip(prctile(y_scen(:,var_plt,:),100-alphau/2,[2 3]))],'--','LineWidth',2,'FaceColor',dflt_col{s},'EdgeColor',dflt_col{s},'LineStyle','--','FaceAlpha',.1,'LineWidth',1)

%Plot credible bands for the 'moderate' policy
y_scen = y_count_s_int;
patch([qe_dt ; flip(qe_dt)],[prctile(y_scen(:,var_plt,:),alphau/2,[2 3]) ; flip(prctile(y_scen(:,var_plt,:),100-alphau/2,[2 3]))],'--','LineWidth',2,'FaceColor',dflt_col{s+1},'EdgeColor',dflt_col{s+1},'LineStyle','--','FaceAlpha',.1,'LineWidth',1)

yline(0)
hold off
axis tight
box on
ylim([-0.6 0.605])
set(gca,'TickLabelInterpreter','latex')

%Define the legend
lh = legend(['Realized', string(scen_lbl)],'Orientation','Horizontal','Location','Northoutside','NumColumns',3,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Indicate recessions
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
fz.Children(1).Children(1).String = fz.Children(1).Children(1).String(1:3);
fz.Children(1).Children(1).FontSize = fz.Children(1).Children(1).FontSize + 2;
fz.Children(1).Children(2).FontSize = fz.Children(1).Children(2).FontSize + 2;
fz.Children(1).Children(2).YAxis(2).TickValues = [];
yyaxis(fz.Children(1).Children(2),'left');
fz.Children(1).Children(2).Children(4).Color = dflt_col{3};
fz.Children(1).Children(2).Children(2).FaceColor = dflt_col{3};
fz.Children(1).Children(2).Children(2).EdgeColor = dflt_col{3};
fz.Children(1).Children(2).Children(5).Color = dflt_col{1};
fz.Children(1).Children(2).Children(3).FaceColor = dflt_col{1};
fz.Children(1).Children(2).Children(3).EdgeColor = dflt_col{1};
fz.Children(1).Children(2).Children(6).Color = 'black';

print(fig_pth + "\MAIN_2",'-dpdf','-vector','-bestfit');

%% ---- ---- Web appendix
%% ---- ---- Figure G.10: Posterior distributions of moments in counterfactual policy experiments

%Time periods to exclude in calculation of moments
t_excl = year(qe_dt) == 2020;

%Compute them
yr_mmnt = [mean(Y_EA(~t_excl,:))' std(Y_EA(~t_excl,:))' skewness(Y_EA(~t_excl,:))' kurtosis(Y_EA(~t_excl,:))'-3];
yk_mmnt = nan(numel(scen_lbl),K,4,pDraw);

for s = 1:numel(scen_lbl)
    y_scen = y_count_med*(s==1) + y_count_s_int*(s==2);

    yk_mmnt(s,:,1,:) = mean(y_scen(~t_excl,:,:),1);
    yk_mmnt(s,:,2,:) = std(y_scen(~t_excl,:,:),1);
    yk_mmnt(s,:,3,:) = skewness(y_scen(~t_excl,:,:),1);
    yk_mmnt(s,:,4,:) = kurtosis(y_scen(~t_excl,:,:),1)-3;
end

%Plot the posterior distributions of the moments
mmnt_lbl = {'Mean','Standard deviation','Skewness','Excess kurtosis'};

%Posterior distribution of moments
figure
tiledlayout(3,K,"TileSpacing","compact"',"Padding","compact")
set(gcf,'Units','centimeters')
figp = get(gcf,'Position');
figp(3) = 29.7;
figp(4) = 21*0.9;
set(gcf,'Position',figp);
set(gcf,'PaperOrientation','landscape')
for mmnt = 1:3
    for k = 1:K
        nexttile
        hold on
        xline(yr_mmnt(k,mmnt),'-','Color',dflt_col{1},'LineWidth',1.5)
        xlm = nan(2,1);
        mdx = [];
        for s = 1:2
            [f,xi,bw] = ksdensity(squeeze(yk_mmnt(s,k,mmnt,:)),'NumPoints',5000);
            patch([xi' ; flip(xi')],[f' ; zeros(size(f'))],'-','LineWidth',2,'FaceColor',dflt_col{s+1},'EdgeColor',dflt_col{s+1},'LineStyle','-','FaceAlpha',.1,'LineWidth',1);
            %Compute percentiles of the smoothed density
            CDF_prb = cumsum(f')/sum(f);
            PDF_prb = f'/sum(f);
            xlm(1) = min(xi(CDF_prb < 0.1 & lagmatrix(PDF_prb,-1) + CDF_prb >= 0.1),xlm(1));
            xlm(2) = max(xi(CDF_prb < 0.9 & lagmatrix(PDF_prb,-1) + CDF_prb >= 0.9),xlm(2));
            mdx = [mdx min(xi(f == max(f)))];
        end
        axis tight
        xlm(1) = min(yr_mmnt(k,mmnt),xlm(1));
        xlm(1) = (1-0.05*sign(xlm(1)))*xlm(1);
        xlm(2) = max(yr_mmnt(k,mmnt),xlm(2));
        xlm(2) = (1+0.05*sign(xlm(2)))*xlm(2);
        xlim(xlm')
        for s = 1:numel(scen_lbl)
            xline(mean(squeeze(yk_mmnt(s,k,mmnt,:))),'--','Color',dflt_col{s+1},'LineWidth',1.5)
        end
        hold off
        box on
        if k == 1
            ylabel(mmnt_lbl{mmnt},'interpreter','latex')
        end
        if mmnt == 1
            title(ttl{k},'Interpreter','latex')
        end
        set(gca,'TickLabelInterpreter','latex')
    end
end

%Define the legend
lh = legend(["Realised",scen_lbl],'Orientation','Horizontal','Location','Northoutside','NumColumns',numel(scen_lbl)+2,'interpreter','latex');
lh.Layout.Tile = 'North';

%Call function to format active figure
PlotFmt([],1)

%Ad hoc modifications
fz = gcf();
fz.Children(1).Children(1).FontSize = fz.Children(1).Children(1).FontSize + 1;
for f = 1:length(fz.Children(1).Children)
    if strcmp(fz.Children(1).Children(f).Type,'axes')
        yyaxis(fz.Children(1).Children(f),'left');
        fz.Children(1).Children(f).Children(2).Color = dflt_col{1};
        fz.Children(1).Children(f).Children(3).FaceColor = dflt_col{3};
        fz.Children(1).Children(f).Children(3).EdgeColor = dflt_col{3};
        fz.Children(1).Children(f).Children(5).Color = 'black';
        fz.Children(1).Children(f).YAxis(1).TickValues = [];
        fz.Children(1).Children(f).YAxis(2).TickValues = [];
        fz.Children(1).Children(f).FontSize = fz.Children(1).Children(f).FontSize + 1;
    end
end

print(fig_pth+"\WA_G10",'-dpdf','-vector','-bestfit');

%% ---- Table 1: Moments of euro area macro variables across policy experiments

%Credible level
alphau = 5;

%Create a corresponding table
sri_tbl = array2table([
    cellstr(string(yr_mmnt(sri_o,:)) + newline + "[]")
    cellstr(string(squeeze(mean(yk_mmnt(:,sri_o,:,:),[2 4]))) + newline + "[" + string(squeeze(prctile(yk_mmnt(:,sri_o,:,:),alphau/2,[2 4]))) + "," + string(squeeze(prctile(yk_mmnt(:,sri_o,:,:),100-alphau/2,[2 4]))) + "]")
    ],'VariableNames',mmnt_lbl);
inf_tbl = array2table([
    cellstr(string(yr_mmnt(inf_o,:)) + newline + "[]")
    cellstr(string(squeeze(mean(yk_mmnt(:,inf_o,:,:),[2 4]))) + newline + "[" + string(squeeze(prctile(yk_mmnt(:,inf_o,:,:),alphau/2,[2 4]))) + "," + string(squeeze(prctile(yk_mmnt(:,inf_o,:,:),100-alphau/2,[2 4]))) + "]")
    ],'VariableNames',mmnt_lbl);
gdp_tbl = array2table([
    cellstr(string(yr_mmnt(gdp_o,:)) + newline + "[]")
    cellstr(string(squeeze(mean(yk_mmnt(:,gdp_o,:,:),[2 4]))) + newline + "[" + string(squeeze(prctile(yk_mmnt(:,gdp_o,:,:),alphau/2,[2 4]))) + "," + string(squeeze(prctile(yk_mmnt(:,gdp_o,:,:),100-alphau/2,[2 4]))) + "]")
    ],'VariableNames',mmnt_lbl);
ciss_tbl = array2table([
    cellstr(string(yr_mmnt(ciss_o,:)) + newline + "[]")
    cellstr(string(squeeze(mean(yk_mmnt(:,ciss_o,:,:),[2 4]))) + newline + "[" + string(squeeze(prctile(yk_mmnt(:,ciss_o,:,:),alphau/2,[2 4]))) + "," + string(squeeze(prctile(yk_mmnt(:,ciss_o,:,:),100-alphau/2,[2 4]))) + "]")
    ],'VariableNames',mmnt_lbl);
irt_tbl = array2table([
    cellstr(string(yr_mmnt(irt_o,:)) + newline + "[]")
    cellstr(string(squeeze(mean(yk_mmnt(:,irt_o,:,:),[2 4]))) + newline + "[" + string(squeeze(prctile(yk_mmnt(:,irt_o,:,:),alphau/2,[2 4]))) + "," + string(squeeze(prctile(yk_mmnt(:,irt_o,:,:),100-alphau/2,[2 4]))) + "]")
    ],'VariableNames',mmnt_lbl);

mmnt_tbl = table(sri_tbl,inf_tbl,gdp_tbl,ciss_tbl,irt_tbl,'RowNames',["Realised" scen_lbl],'VariableNames',ttl);

%% END