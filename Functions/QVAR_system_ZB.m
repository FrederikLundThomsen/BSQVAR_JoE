% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

function [D_p,A0_p,AL_p,BL_p,V_p,eps_p,lambdaT,sigmaT,nu_p] = QVAR_system_ZB(Y,P,X,S,D,Quant,ZeroRestrictions,niter,Bprior,VBprior,tightness)

%Estimate and store the parameters at the grid of quantiles
[T,K] = size(Y);

[~,C] = size(D);

[~,n] = size(Quant);

if isempty(X) == 0
    M = size(X,2);
else
    M = 0;
end

%Pre-allocate memory to store coefficients
D_p = zeros(K,C,n,niter);
A0_p = zeros(K,K,n,niter);
if P > 0
    AL_p = zeros(K,K*P,n,niter);
else
    AL_p = [];
end
if ismissing(X) == 0
    BL_p = zeros(K,M*(S+1),n,niter);
else
    BL_p = [];
end

%Posterior covariance matrix of parameters
% K endogenous varialbes x C+K*P+M*(S+1)+K parameters x # of quantiles (n) x posterior draws (niter)
V_p = zeros(K,C+K*P+M*(S+1)+K,C+K*P+M*(S+1)+K,n,niter);

%Initialise tightness / looseness paramter. One for each variable at each
%quantile
if tightness == 1
    lambdaT = zeros(K,n,niter);
else
    lambdaT = [];
end

%Array to store latent variables
nu_p = zeros(K,T-P,n,niter);

%Prepare estimation sample
if isempty(ZeroRestrictions) == 1
    ZeroRestrictions = ones(C+K*P+M*(S+1)+K,K);
end

%Z is the matrix of regressors
%Feed it deterministic terms
Z = D(P+1:end,:);

%Endogenous lags (Simultaneous regressors are added later)
if P > 0
    for p=1:P
        Z = [Z Y((P-p+1):(end-p),:)];
    end
end

%Add exogenous lags and simultaneous regressors if necessary
if isempty(X) == 0    
    for s=0:min(P,S)
        if s == 0
            Xl = X((P+1):end,:);
        else
            Xl = [Xl X((P+1-s):(end-s),:)];
        end
        if s == min(P,S)
            Z = [Z Xl];
        end
    end
end

%Array to store errors
eps_p = zeros(size(Z,1),K,n,niter);

%Array to store volatility
sigmaT = zeros(K,T-P,n,niter);

%Extract the priors
Bprior_e = cell2mat(Bprior);
VBprior_e = cell2mat(VBprior);

%Initialise the relevant starting values
nu = ones(T-P,1);
sigma2 = zeros(K,T-P,n);

tau_sq = 2./(Quant.*(1-Quant));
theta = (1-2*Quant)./(Quant.*(1-Quant));

%Hyper parameters for scaling parameter lambda (Prior mean is 3, mode is 1.5 and
%std is 3)
alpha_lamdba = 3;
gamma_lambda = 6;

%Hyperparameters on the inverse-gamma prior for the standard error of the error term
alpha_sigma = 0.01;
gamma_sigma = 0.01;

%Start a timer
TIME = tic;

%Create cell with estimating data. 
% We use cells because the number of
%regressors differs for each endogenous variable depending on the ordering
%as well as additional zero restrictions
Zr = {};
Zempty = {};
for i = 1:K
    %Temporary vector with zero restrictions for variable i
    Zeroi = ZeroRestrictions(:,i);
    %For the first endogenous variable, with no simultaneous endogenous
    %variables on the right hand side
    if i == 1
        %Set restricted columns to zero
        Zir = Z.*Zeroi(1:C+K*P+M*(S+1))';
        %Store the location of the Z columns that are removed
        Zempty = [Zempty {find(~any(Zir,1))}];
        %Remove any zero columns/ regressors
        Zir(:,~any(Zir,1)) = [];    
    %For the remaining endogenous variables
    else
        %Add the relevant simultaneous observations of endogenous variables
        %(1:(i-1))
        Zir = Z;
        for j=1:(i-1)
            Zir=[Zir,Y((P+1):end,j)];
        end
        
        %Remove restricted columns
        Zir = Zir.*Zeroi(1:C+K*P+M*(S+1)+(i-1))';
        Zempty = [Zempty {find(~any(Zir,1))}];
        Zir(:,~any(Zir,1)) = [];    
    end
    Zr = [Zr Zir];
end

warning('off')

%Estimate system equation by equation for each MC draw
for i=1:K
        
    %Retrieve the set of regressors
    Ze = Zr{i};
    %And the location of zero restricted regressors - To be used to add
    %zeros in the posterior draw of parameters and its covariance matrix
    Zempte = Zempty{i};
    
    %Loop over quantiles
    for j=1:n

        %Select the prior parameters relevant for the given
        %variable/quantile

        %Simultaenous endogenous observations are ordered last, hence we
        %can leave out the last K-(i-1) prior means and variances
        Bij_prior = squeeze(Bprior_e(i,1:end-(K-(i-1)),j))';
        VBij_prior = squeeze(VBprior_e(i,1:end-(K-(i-1)),j))';

        %Impose zero restrictions by removing restricted regressors
        Bij_prior(Zempte) = [];
        VBij_prior(Zempte) = [];
        
        %Initialise parameters
        bz_ij = Bij_prior;

        s_b = length(bz_ij);

        %Prior mean and variance of coefficients in the univariate regression
        mu_beta = Bij_prior;
        Sigma_beta_inv = inv(diag(VBij_prior));
        
        %Initialize time series of latent variable
        nu_ij = nu;

        %Initialise looseness.
        lambda = 1e+0;
        
        %Start the Gibbs sampler
        for ii = 1:niter
        
            %Print status at every 500th posterior draw
            if mod(ii,500) == 0
                "[" + string(datetime()) + "] " + "Time elapsed: " + sprintf("%s",duration([0, 0, toc(TIME)])) + ". Variable: " + i + " Quantile: "  + floor(round(Quant(j),2)*100) + " Iteration: " + ii + "."
            end
            
            %Auxiliary variable
            y_tilde = Y(P+1:end,i) - theta(j)*nu_ij;

            % Sample variance from its posterior
            alpha_sigma_bar = alpha_sigma + 3*(T-P)/2;
            gamma_sigma_bar = gamma_sigma + sum(nu_ij) + 0.5*sum((y_tilde - Ze*bz_ij).^2./(tau_sq(j)*nu_ij));       

	        %Draw from the IG(alpha_sigma_bar,gamma_sigma_bar) distribution, using that if X ~ G(a,1) then Y = 1/X ~ IG(a,1), while if Y ~ IG(a,b) then kY ~ IG(a,kb).
            sigma2(i,:,j) = gamma_sigma_bar/randg(alpha_sigma_bar);

            %Compute posterior covariance matrix of regression coefficients
            U = diag(1./(sqrt(sigma2(i,:,j)'.*tau_sq(j).*nu_ij)));
            Z_tilde = (Ze'*U)';
            Y_tilde = (y_tilde'*U)';
            Sigma_beta_bar = inv(Z_tilde'*Z_tilde + Sigma_beta_inv./lambda);

            %Determine posterior mean of regression coefficients
            mu_beta_bar = Sigma_beta_bar*(Z_tilde'*Y_tilde + Sigma_beta_inv./lambda*mu_beta);    

            %Draw regression coefficients from their posterior
            bz_ij = mu_beta_bar + chol(Sigma_beta_bar)'*randn(s_b,1);
        
            % Sample latent variables nu(t)
            for t = 1:T-P            
                k1 = sqrt(theta(j).^2 + 2*tau_sq(j))/abs(Y(t+P,i)-Ze(t,:)*bz_ij);
                k2 = (theta(j).^2 + 2*tau_sq(j))./(sigma2(i,t,j)*tau_sq(j));
                nu_ij(t) = 1./igrnd(k1,k2); 
            end
            
            if tightness == 1
                %Draw the EA scaling parameter from its posterior
                alpha_lambda_bar = alpha_lamdba+s_b/2;
                gamma_lambda_bar = (gamma_lambda + 0.5*(bz_ij - mu_beta)'*Sigma_beta_inv*(bz_ij - mu_beta));
                lambda = gamma_lambda_bar/randg(alpha_lambda_bar);

                lambdaT(i,j,ii) = lambda;
            end

            %Reinstate the relevant zeros in the complete vector of
            %parameter draws, i.e. including restricted regressors
            bi = zeros(C+K*P+M*(S+1)+i-1,1);
            bi(~ismember((1:C+K*P+M*(S+1)+i-1),Zempte)) = bz_ij;

            %Same for the posterior covariance matrix of parameters
            Vi = zeros(C+K*P+M*(S+1)+K,C+K*P+M*(S+1)+K);
            Vi(~ismember(1:C+K*P+M*(S+1)+i-1,Zempte),~ismember(1:C+K*P+M*(S+1)+i-1,Zempte)) = Sigma_beta_bar;

            %Store the posterior parameter draws in the designated arrays

            %Deterministic terms
            D_p(i,:,j,ii) = bi(1:C,:);

            %Simultaneous effects
            if i > 1
                A0_p(i,1:i-1,j,ii) = bi(end-(i-2):end,:);
            end
            
            %Endogenous lags
            if P > 0
                %Endogenous lags
                AL_p(i,:,j,ii) = bi(1+C:(C+K*P),:);
            end
        
            %Exogenous variables
            if isempty(X) == 0
                BL_p(i,:,j,ii) = bi(1+C+K*P:C+K*P+M*(S+1),:);
            end

            %Store other relevant posterior outputs
            V_p(i,:,:,j,ii) = Vi;

            eps_p(:,i,j,ii) = Y(P+1:end,i) - Ze*bz_ij;
 
            sigmaT(i,:,j,ii) = sigma2(i,:,j);

            nu_p(i,:,j,ii) = nu_ij;

        end
        
    end
   
end
warning('on')

end