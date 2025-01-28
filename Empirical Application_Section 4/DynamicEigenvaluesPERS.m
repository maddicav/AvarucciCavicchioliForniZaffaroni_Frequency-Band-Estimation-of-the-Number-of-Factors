%PERIODOGRAM SMOOTHING ESTIMATOR
%
function [lambda, lambdam]  = DynamicEigenvaluesPERS(X,neigs,M)

[T,~] = size(X);

if nargin < 3
    M = round(.75*sqrt(T));
end
if nargin < 2
    neigs = 2*M+1;
end

W = 2*M + 1;
lambda = zeros(floor((T-1)/2),W);
opts.disp=0;
% compute the eigenvalues  for floor((T-1)/2)+1 points 2*pi*h/T, h=0,...,floor((T-1)/2)
for h = 0:floor((T-1)/2)
    approx = h-M:h+M;
Xf=X'*exp(-sqrt(-1)*(1:T)'*2*pi*approx/T);
lambda(h+1,:) = abs(eigs(Xf*Xf'/(W*T),W,'lr',opts));

end
lambdam = lambda(:,end);
lambda = lambda(:,1:neigs);