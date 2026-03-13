function y = igrnd(mu,lambda)
% Generate a draw from the Inverse Gaussian (Wald) distribution
%----------------------------------------------------------------
% This method is from
% John R. Michael, William R. Schucany, and Roy W. Haas (1976) Generating 
% Random Variates Using Transformations with Multiple Roots. The American 
% Statistician 30: 88-90.
%

T = length(mu);

% First generate a chi_square(1) variate
v0 = randn(T,1).^2;

% intermediate steps
w = mu.*v0;
c = mu./(2*lambda);

% Now find the two roots in Eq. (5) of the paper
x1 = mu+c.*(w-sqrt(w.*(4*lambda+w)));
x2 = (mu.^2)./x1;

% Pick x2 if Bernoulli trial fails
p1_v0 = mu./(mu+x1);
y = x1;
idx = rand(T,1)>p1_v0;
y(idx) = x2(idx);

end
