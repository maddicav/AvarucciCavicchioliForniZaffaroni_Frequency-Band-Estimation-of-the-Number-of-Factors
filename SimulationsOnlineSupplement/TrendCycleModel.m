%--------------------------------------------------------------------------
% variables = TrendCycleModel(N,T,q,s)
% INPUT    N,T    :  panel dimension
%          q      :  number of factors/shocks
%          s      : standard deviation idiosyncratic component
% OUTPUT   data:  : (T x n) panel data (see fourth experiment, Trend-Cycle Model)
%--------------------------------------------------------------------------
%
function variables = TrendCycleModel(N,T,q,s)
if nargin <4
    s=1;
end
%
u = randn(T+50,q+1);
%--------------------------------------------------------------------------
% Coefficients Filters Permanent Shock
%--------------------------------------------------------------------------
coeffperm0 = rand(N,q)*2-1;
coeffperm1 = rand(N,q)*1.6-1;
%--------------------------------------------------------------------------
% Coefficients Filters Transitory Shock
%--------------------------------------------------------------------------
coefftrans1 = rand(N,1)*2-1;
coefftrans2 = rand(N,1)*.7;
%
common=zeros(T,N); % Preallocation
%
for i=1:N
   a = zeros(T+50,N);
   
      for j=1:q
         
         a = a + filter(coeffperm0(i,j),[1  -coeffperm1(i,j)],u(:,j));
         
      end
      a = a + filter(coefftrans1(i)*[1  -1],[1 -coefftrans2(i)],u(:,q+1));
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



