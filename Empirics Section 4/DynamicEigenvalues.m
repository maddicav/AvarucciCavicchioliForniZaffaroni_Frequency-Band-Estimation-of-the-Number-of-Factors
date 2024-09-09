%PERIODOGRAM SMOOTHING ESTIMATOR
%
function [D,V, Sigma]  = DynamicEigenvalues(X,q,M)

[T,n] = size(X);


if nargin < 3
    M = round(.75*sqrt(T));
end
% if nargin < 2
%     neigs = 2*M+1;
% end
W = 2*M + 1;
if nargin < 2
q = W;
end
nof = floor((T-1)/2)+1;


% compute the spectral density matrix and the eigenvalues  for floor((T-1)/2)+1 
%points 2*pi*h/T, h=0,...,floor((T-1)/2)
Sigma = zeros(n,n,nof);
V = zeros(n,q,nof);
D = zeros(q,q,nof);
for h = 0:nof-1
    approx = h-M:h+M;
Xf=X'*exp(-sqrt(-1)*(1:T)'*2*pi*approx/T);
Sigma(:,:,h+1) = Xf*Xf'/(W*T);
[V(:,:,h+1),D(:,:,h+1)]  = eigs(Sigma(:,:,h+1),q,'lr');


end
