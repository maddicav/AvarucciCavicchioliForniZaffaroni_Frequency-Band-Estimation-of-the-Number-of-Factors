% Compute the eigevalues of the smoothed periodogram for a sequence X
%
% [lambda, lambdam]  = DynamicEigenvaluesPERS(X,neigs,M)
% INPUT:  X         : T x n data matrix  (required)
%         neigs     : upper bound on the number of factors (default 2M+1)  
%         M         : the bandwidth is computeas as M=[c(sqrt(T)) ] (default M=[.75*sqrt(T)])
%
% OUTPUT: lambdam   : smallest eigenvalue
%         lambda    : neigs eigenvalues, dereasing order
%--------------------------------------------------------------------------
function [lambda, lambdam]  = DynamicEigenvaluesPERS(X,neigs,M)
%
[T,n] = size(X);
%
if nargin < 3
    M = round(.75*sqrt(T));
end
if nargin < 2
    neigs = 2*M+1;
end
%
W = 2*M + 1;
lambda = zeros(floor((T-1)/2),min(n,W));
opts.disp=0;
for h = 0:floor((T-1)/2)
    approx = h-M:h+M;
Xf=X'*exp(-sqrt(-1)*(1:T)'*2*pi*approx/T);
lambda(h+1,:) = abs(eigs(Xf*Xf'/(W*T),min(n,W),'lr',opts));
end
lambdam = lambda(:,end);
lambda = lambda(:,1:neigs);