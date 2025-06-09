% log criterion to determine the number of dynamic factors according to Hallin and Liska (2007) 
% computed with penalty p_1 and M_T=[.75*sqrt(T)]
%
% qlog=HLcrit(panel, qmax,cmax, ngrid,Tgrid) 
% INPUT:    panel           :   T x n data matrix  (required)
%           qmax            :   upper bound on the number of factors (required)  
%           ngrid, Tgrid    :   row vectors of dimension "J" with elements
%                               n_1<n_2<...<n_J=n and T_1<T_2<...<T_J=T
%                               T_j x n_j subpanels are used where
%                               (default values: ngrid = n-30:10:n; Tgrid =T-30:10:T.                            
%           cmax           :   c = [0:cmax] (default value: 3)
%     
% OUTPUT:   qlog           :  estimated number of factors
%                               
%           
%----------------------------------------------------------------------------------------------
function qlog = HLcrit(panel, qmax,cmax, ngrid,Tgrid)
% 
[T,n] = size(panel);
%
if nargin <= 2    
    ngrid = n-30:10:n;
    Tgrid = T-30:10:T;
    cmax = 3;
end
% Preallocation
q=zeros(length(ngrid),cmax*100);
%
for j=1:length(ngrid)
    
    
    TT = Tgrid(j);
    nn = ngrid(j);
    window = floor(sqrt(TT)*0.75);
    eigv = subr_1(panel(1:TT,1:nn), window); % see the function E = subr_1(x, M) below
    IC1 = flipud(cumsum(flipud(eigv)));
    IC1 = IC1(1:qmax+1,:);

    p = ((window/TT)^0.5 + window^(-2) + nn^(-1))*log(min([(TT/window)^0.5;  window^2; nn]));
    
    
    for c = 1:cmax*100
        cc = c/100;
       IC_log = log(IC1/nn) + (0:qmax)'*p*cc;
       [~,qq] = min(IC_log);
       q(j,c) = qq-1;
        
       
    end %c
end %n


qstableindex = find(std(q) == 0);
qstableindexnoqmax = find(q(end, qstableindex) < qmax);
qlog = q(end, qstableindex(qstableindexnoqmax(1)) );
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% 
%This function compute of the eigenvalue of the lag window estimator of the
%spectral density. See equation (6) in Hallin and Liska (2007)
%------------------------------------------------------------------------
function  E = subr_1(x, M)
[T,N] = size(x);
w = M;
W = 2*w+1;
% compute covariances 
SS = zeros(W,N*N);
for k = 1:w+1
     Sk = center(x(k:T,:))'*center(x(1:T+1-k,:))/(T-k);
     S_k = Sk';
     SS(w-k+2,:) = S_k(:)';  
     SS(w+k,:) = Sk(:)';
end


% compute the spectral matrix in w points (S)


freq = 0:2*pi/W:pi;
Factor = exp(-sqrt(-1)*(-w:w)'*freq);

ww = 1 - abs(-w:w)/(w+1);

% multiply first line of SS by ww(1)
% multiply second line of SS by ww(2)

S = diag(ww)*SS(1:W,:); 

% take the inverse of the vec operator
% S = reshape(S'*Factor,N,N*W);

S = reshape(S'*Factor,N,N*(w+1));

% compute the eigenvalues 
eigenvalues = zeros(N,w+1);
D = eig(S(:,1:N));
eigenvalues(:,1) = sort(real(D),'descend');

for j = 1:w
   D = eig(S(:,j*N+1:(j+1)*N));
   eigenvalues(:,1+j) = sort(real(D),'descend');
end
%eigenvalues
%[eigenvalues(:,1)  eigenvalues(:,2:jj+1)*2]
E = [eigenvalues(:,1)  eigenvalues(:,2:w+1)*2]*ones(w+1,1)/(2*w+1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XC = center(X)

[T , ~] = size(X);
XC = X - ones(T,1)*mean(X); 
XC = XC./kron(ones(T,1),std(XC));





