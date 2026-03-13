% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

function [Y_F,X_F] = QVAR_Forecast_Scenario_Z_FixP(Y,P,X,S,D,Quant,fDt,H,N,ScenarioDesign,ZeroRestrictions,DC,A0,AL,BL,BX,seed)

%Monte Carlo forecasts of the QVAR using fixed parameters

%Initialize relevant parameters
[~,K] = size(Y);

[~,n] = size(Quant);

C = size(D,2);

if isempty(X) == 0
    M = size(X,2);
else
    M = 0;
end

%Prepare estimation sample
if isempty(ZeroRestrictions) == 1
    ZeroRestrictions = ones(C+K*P+M*(S+1)+K,K);
end

%Compute weighted probabilities if the set of quantiles is not equidistant
wq = [];
for q = 1:n
    if n == 1
        wq = 1;
    else
        if q == 1
            wq = [wq; Quant(q) + (Quant(q+1) - Quant(q))/2];
        elseif q > 1 & q < n
            wq = [wq; (Quant(q+1) - Quant(q-1))/2];
        else 
            wq = [wq;1-Quant(q) + (Quant(q) - Quant(q-1))/2];
        end
    end
end

%Apply the seed if relevant
if ~isempty(seed)
    rng(seed);
end

%Draw all random quantiles to be used for the forecasting
rQuants = randsample(1:n,(K+M)*H*N,true,wq);
rQuants = reshape(rQuants,K+M,H,N);

%Number of quantile decimals
qdec = cellfun(@(c)strsplit(c,'.'),string(Quant),'UniformOutput',false);
qdeca = 0;
for d = 1:length(qdec)
    if length(qdec{d}) == 2
        qdeca = max(qdeca,strlength(qdec{d}(2)));
    end
end

%Extract scenarios
if isempty(ScenarioDesign{1}) == 1
    ImpulseDesign = zeros(H,K+M,n);
else
    ImpulseDesign = ScenarioDesign{1};
end
if isempty(ScenarioDesign{2}) == 1
    QuantileDesign = nan(H,K+M);
else
    QuantileDesign = ScenarioDesign{2};
end
if isempty(ScenarioDesign{3}) == 1
    PathDesign = nan(H,K+M,n);
else
    PathDesign = ScenarioDesign{3};
end

%Note that path design sets variables at fixed values and will always
%override Impulse and Quantile designs. Impulse and quantile designs can
%activate at the same time, however.

%Initialize array to hold simulation output and parameter estimates
Y_F = cell(N,1);

if isempty(X) == 0
    X_F = cell(N,1);
else
    X_F = [];
end

%Set restricted coefficients equal to zero
for j = 1:K
    DC(j,ZeroRestrictions(1:C,j) == 0,:) = 0;
    A0(j,ZeroRestrictions(end-(K-1):K,j) == 0,:) = 0;
    AL(j,ZeroRestrictions(C+1:C+K*P,j) == 0,:) = 0;
    if isempty(X) == 0
        BL(j,ZeroRestrictions(C+K*P+1:C+K*P+M*(S+1),j) == 0,:) = 0;
    end
end

%Simulate paths up to fDt+H
parfor i=1:N
    warning('off','all')
    
    %Prints every 1000th iteration to the command window
    if mod(i,1e+3) == 0
        % sprintf('Simulation: %d', i)
    end
    
    Y_N = zeros(H+1,K);
    X_N = zeros(H+1,M);
    
    %Store the last actual observations of Y and X
    Y_N(1,:) = Y(fDt,:);
    
    if isempty(X) == 0
        X_N(1,:) = X(fDt,:); 
    else
        X_N = 0;
    end

    %Create vectors with lagged values of endogenous and exogenous variables
    YL_i = zeros(K*P,1);
    
    %Enrich them from observed data
    for p=1:P
        YL_i((1+K*(p-1)):(K*p)) = Y(fDt-(p-1),:);
    end
    
    if isempty(X) == 0
        XL_i = zeros(M*(P+1),1);
        
        for s=0:P
            XL_i((1+M*s):(M*(s+1))) = X(fDt-s,:);
        end
    else
        XL_i = 0;
    end

    for h=1:H
        %Pick random quantiles from the quantile space
        r_quants = Quant(rQuants(1:K,h,i));
        
        %Deterministic terms
        D_i = zeros(K,C);
        
        %Simultaneous effects
        A0_i = zeros(K);

        %Lag coefficients for Y
        AL_i = zeros(K,K*P);

        %Coefficients for exogenous variables
        if isempty(X) == 0
            BL_i = zeros(K,M*(S+1));
        else
            BL_i = 0;
        end

        %If the forecasting is in-sample, the actual series of deterministic terms is used
        Dt = [1,zeros(1,C-1)];

        %Vector to store predetermined path at horizon h
        Ypath = zeros(K,1);

        %Vector to store quantile shocks
        shk = [];

        %Construct companion form
        for j=1:K
            
            %If the design matrix contains a non-NaN value for variable j
            %for horizon h, then the randomly assigned quantile is replaced
            %with the quantile value specified in the design matrix
            if (sum(~isnan(QuantileDesign(:,j))) > 0)&&(ismember(h,find(~isnan(QuantileDesign(:,j)))))
                r_quants(j) = QuantileDesign(h,j);
            end
            
            %Index for the randomly picked quantile
            qidx = round(r_quants(j),qdeca) == round(Quant,qdeca);
            
            %Insert coefficients in A0
            A0_i(j,:) = squeeze(A0(j,:,qidx));

            %And the rest          
            AL_i(j,:) = squeeze(AL(j,:,qidx));

            if isempty(X) == 0
                BL_i(j,:) = squeeze(BL(j,:,qidx));
            else
                BL_i = 0;
            end

            D_i(j,:) = DC(j,:,qidx);

            %Identify the quantile shocks for the given quantile draw for
            %variable j
            shk = [shk;squeeze(ImpulseDesign(h,j,qidx))]; 
            
            %If the j'th variable is preset at the given horizon, all
            %related rows are set to zero, as it is unaffected by
            %endogenous as well as exogenous variables
            if (sum(~isnan(PathDesign(:,j,qidx))) > 0)&&(ismember(h,find(~isnan(PathDesign(:,j,qidx)))))
                A0_i(j,:) = 0;
                AL_i(j,:) = 0;
                if isempty(X) == 0
                    BL_i(j,:) = 0;
                else
                    BL_i = 0;
                end
                D_i(j,:) = 0;
                
                Ypath(j) = PathDesign(h,j,qidx);
            end
        end
        
        %Forecast X for fDt+h
        if isempty(X) == 0
            %Pick random quantile to forecast X
            x_quants = Quant(rQuants(K+1:K+M,h,i));
            
            %Auxiliary variable
            stp = 0;

            %Initialize vector to store forecast of X
            Xf = [];
            
            %Forecast each individual exogenous variable
            for x=1:M                    
                %Replace random quantil with predetermined quantile if
                %relevant
                if (sum(~isnan(QuantileDesign(:,K+x))) > 0)&&(ismember(h,find(~isnan(QuantileDesign(:,K+x)))))
                    x_quants(x) = QuantileDesign(h,K+x);
                end
                
                %Number of coefficients for each quantile regression
                Ncx = C+P;
                
                %Auxiliary variable, cf. above
                idx = stp+(find(round(Quant,qdeca) == round(x_quants(x),qdeca))-1)*Ncx;                
                
                %Forecast
                Xfm = Dt*BX(idx+1:idx+C) + BX(idx+C+1:idx+C+P)'*XL_i(x:M:end-M) + ImpulseDesign(h,K+x,find(round(Quant,qdeca) == round(x_quants(x),qdeca)));
 
                %If X(m) is predetermined by an exogenous path at the given
                %horizon, this is inserted instead
                if (sum(~isnan(PathDesign(:,K+x,find(round(Quant,qdeca) == round(x_quants(x),qdeca))))) > 0)&&(ismember(h,find(~isnan(PathDesign(:,K+x,find(round(Quant,qdeca) == round(x_quants(x),qdeca)))))))
                    Xfm = PathDesign(h,K+x,find(round(Quant,qdeca) == round(x_quants(x),qdeca)));                    
                end
                
                Xf = [Xf ; Xfm];
                
                %Update auxiliary variable
                stp = stp + n*Ncx;                
            end 

            %Update the vector with exogenous variables and their lags
            XL_i = [Xf ; XL_i(1:end-M)];
            
            %Store forecast of X
            X_N(h+1,:) = Xf;
        end
        
        %Forecast Y for fDt+h and store the output
        %Only the constants from the deterministic terms are picked (Needs
        %to be updated and generalized if a deterministic trend is included)
        if isempty(X) == 0
            Y_N(h+1,:) = (eye(K)-A0_i)\(D_i*Dt' + AL_i*YL_i + BL_i*XL_i(1:(S+1)*M) + Ypath + shk);
        else
            Y_N(h+1,:) = (eye(K)-A0_i)\(D_i*Dt' + AL_i*YL_i + Ypath + shk);
        end
        
        %Update lagged values for endogenous and exogenous variables
        if P == 1
            %If there's only one lag, then this is replaced by the forecast
            YL_i = Y_N(h+1,:)'; 
        else
            %Otherwise the forecast enters the vector of lags while the
            %oldest values are discarded
            YL_i = [Y_N(h+1,:)'; YL_i(1:end-K)]; 
        end
    end
    
    Y_F{i} = Y_N;
    if isempty(X) == 0
        X_F{i} = X_N;
    end
    
end
