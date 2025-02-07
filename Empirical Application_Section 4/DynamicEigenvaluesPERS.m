% Compute the eigevalues of the smoothed periodogram for a sequence X
%
% [lambda, lambdam]  = DynamicEigenvaluesPERS(X,neigs,M)
% INPUT:  X         : T x n data matrix  (required)
%         neigs     : number of largest eignevalues to be returned  (default 2M+1)  
%         M         : the bandwidth is computeas as M=[c(sqrt(T)) ] (default M=[.75*sqrt(T)])
%
% OUTPUT: lambdam   : column vector, j-th entry is the "neigs" eigenvalue 
%                     of the smoothed periodogram at the j-th Fourier Frequency
%                     
%               
%         lambda    : matrix, (j,k) entry is the k (k=1,...,negis) largest eigenvalue 
%                     of the smoothed periodogram at the j-th Fourier Frequency
%--------------------------------------------------------------------------
function [lambda, lambdam]  = DynamicEigenvaluesPERS(X,neigs,M)

[T,~] = size(X);

if nargin < 3
    M = round(.75*sqrt(T));
end
if nargin < 2
    neigs = 2*M+1;
end

W = 2*M + 1;
lambda = zeros(floor((T-1)/2)+1,W);
opts.disp=0;
% compute the eigenvalues  for floor((T-1)/2)+1 points 2*pi*h/T, h=0,...,floor((T-1)/2)
for h = 0:floor((T-1)/2)
    approx = h-M:h+M;
Xf=X'*exp(-sqrt(-1)*(1:T)'*2*pi*approx/T);
lambda(h+1,:) = abs(eigs(Xf*Xf'/(W*T),W,'lr',opts));

end
lambdam = lambda(:,end);
lambda = lambda(:,1:neigs);