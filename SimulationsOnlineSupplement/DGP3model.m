%--------------------------------------------------------------------------
% data = DGP3model(N,T,q,opt) 
% INPUT    N,T    :  panel dimension
%          q      :  number of factors/shocks
%          opt    :  'LI'or 'SI'(Large or small idosyncratic component)
% OUTPUT   data:  :  (T x n) panel data (see third experiment)
%--------------------------------------------------------------------------
function data = DGP3model(N,T,q,opt)
%--------------------------------------------------------------------------
% Idiosyncratic component
%--------------------------------------------------------------------------
v = zeros(T+50,N);
e = v;
uu = randn(T+50,N);
v(:,1) = uu(:,1);
for i = 2:N
    v(:,i) = 0.2*v(:,i-1) + uu(:,i);
end
rho = rand(1,N)*1.6-0.8;
e(1,:) = v(1,:);
for t = 2:T+50
    e(t,:) = rho.*e(t-1,:) + v(t,:);
end
idio = e(51:end,:);
%--------------------------------------------------------------------------
% Common component
%--------------------------------------------------------------------------
u = randn(T+50,q);
coeffMA0=rand(q,N)*2-1;
coeffMA1=rand(q,N)*2-1;
coeffMA2=rand(q,N)*2-1;
coeffAR1=rand(q,N)*1.6-.8;
coeffAR2=rand(q,N)*1.6-.8;
common=zeros(T,N);
%
for i=1:N
   a = zeros(T+50,1);  
      for j=1:q
         ARfilt = fliplr(conv([coeffAR1(j,i) 1],[coeffAR2(j,i) 1]));
         a = filter([ coeffMA0(j,i) coeffMA1(j,i) coeffMA2(j,i)],ARfilt,u(:,j)) + a;    
      end
   common(:,i) = a(51:T+50);
end
stdcommon = sqrt(sum(var(common))/N);
common = common/stdcommon;
stdidio = sqrt(sum(var(idio))/N);
idio = idio/stdidio;
if  strcmp(opt,'LI')==1  % large idio
idio = idio*1;
else
    idio = idio*.5;      %small idio
end
data = common+idio;



