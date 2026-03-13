% ---------------------------------------------------------------- %
% Part of the BSQVAR Toolbox for Chavleishvili, S., Engle, R., 
% Fahr, S., Kremer, M. Lund-Thomsen, F., Manganelli, S. and 
% Schwaab, B. (2026). 'Macro-prudential policy under asymmetric risks:
% A Bayesian Structural Quantile VAR approach'. Journal of Econometrics.
%
% Available from GitHub: https://github.com/FrederikLundThomsen/BSQVAR_JoE
% ---------------------------------------------------------------- %

function [b,q,eps,Zr] = VAR_Estim_Z(Y,P,X,S,D,i,ZeroRestrictions)

[T,K]= size(Y);

[~,C] = size(D);

if isempty(X) == 0
    M = size(X,2);
else
    M = 0;
end

if isempty(ZeroRestrictions) == 1
    ZeroRestrictions = ones(C+K*P+M*(S+1)+K*(K-1)/2,K);
end

Zeroi = ZeroRestrictions(:,i);      

Z = D(P+1:end,:);

if P > 0
    for p=1:P
        Z = [Z Y((P-p+1):(end-p),:)];
    end
end

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

if i==1
    %Impose zero restrictions
    Zr = Z.*Zeroi(1:C+K*P+M*(S+1))';
    %Store the location of the Z columns that are removed
    Zempty = find(~any(Zr,1));
    %Remove any potential zero columns
    Zr(:,~any(Zr,1)) = [];

    bz = (Zr'*Zr)\(Zr'*Y((P+1):end,i));

    %Reinstate zeros in b
    b = zeros(size(Z,2),1);
    b(~ismember((1:size(Z,2)),Zempty)) = bz;
    
else
    for j=1:(i-1)
        Z=[Z,Y((P+1):end,j)];
    end
    
    Zr = Z.*Zeroi(1:C+K*P+M*(S+1)+(i-1))';
    Zempty = find(~any(Zr,1));
    Zr(:,~any(Zr,1)) = [];
    
    bz = (Zr'*Zr)\(Zr'*Y((P+1):end,i));
    
    b = zeros(size(Z,2),1);
    b(~ismember((1:size(Z,2)),Zempty)) = bz;
end

q=Z*b;
eps=Y((P+1):end,i)-q;

end
