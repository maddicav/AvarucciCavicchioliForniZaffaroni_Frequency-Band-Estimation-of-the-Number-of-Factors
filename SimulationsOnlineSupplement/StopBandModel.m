%--------------------------------------------------------------------------
% variables = StopBandModel(N,T,q,s)
% INPUT    N,T    :  panel dimension
%          q      :  number of factors/shocks
%          s      : standard deviation idiosyncratic component
% OUTPUT   data:  : (T x n) panel data (see fourth experiment, Stop-Band Model)
%--------------------------------------------------------------------------

function variables = StopBandModel(N,T,q,s)
if nargin <4
    s=1;
end

u = randn(T+50,q+1);
common=zeros(T,N); % Preallocation
%--------------------------------------------------------------------------
% Coefficients Filters First Shock
%--------------------------------------------------------------------------
coeffperm0 = rand(N,q)*2-1;
coeffperm1 = rand(N,q)*1.6-1;
%--------------------------------------------------------------------------
% Coefficients Filters Second Shock
%--------------------------------------------------------------------------
b=conv([1 -exp(1i*2*pi*20/240)],[1 -exp(-1i*2*pi*20/240)]);
coeffperm0(:,q+1) = rand(N,1)-.5;
coeffperm1(:,q+1) = rand(N,1)*.1+.8;
%--------------------------------------------------------------------------
for i=1:N
   a = zeros(T+50,N);
   
      for j=1:q
         
         a = a + filter(coeffperm0(i,j),[1  -coeffperm1(i,j)],u(:,j));
         
      end
      a = a + coeffperm0(i,q+1)*filter(b,[1 -coeffperm1(i,q+1)],u(:,q+1));
   common(:,i) = a(51:T+50);
end
%--------------------------------------------------------------------------
coeffidio = (rand(N,1)*2-1);
idio = randn(T,N)*diag(coeffidio);
stdcommon = sqrt(sum(var(common))/N);
common = common/stdcommon;
stdidio = sqrt(sum(var(idio))/N);
idio = s*idio/stdidio;
%
variables = standardize(common+idio);







