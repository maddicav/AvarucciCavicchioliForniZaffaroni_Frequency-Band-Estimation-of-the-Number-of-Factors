% generate data from
% ar factor model with q factors
%
function variables = modelloARMA(N,T,q,opt1,opt2)
if nargin == 3
    opt1 = 'LI';
    opt2 = 'ARMA';
end
u = randn(T+50,q)*diag(rand(1,q)*0+1);
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
coeffMA0=rand(q,N)*2-1;
coeffMA1=rand(q,N)*2-1;
coeffMA2=rand(q,N)*2-1;
coeffAR1=rand(q,N)*1.6-.8;
coeffAR2=rand(q,N)*1.6-.8;
if strcmp(opt2, 'ARMA')
for i=1:N
   a = zeros(T+50,1);
   
      for j=1:q
         ARfilt = fliplr(conv([coeffAR1(j,i) 1],[coeffAR2(j,i) 1]));
         a = filter([ coeffMA0(j,i) coeffMA1(j,i) coeffMA2(j,i)],ARfilt,u(:,j)) + a;
         
      end
   common(:,i) = a(51:T+50);
end
elseif strcmp(opt2, 'MA')
for i=1:N
   a = zeros(T+50,1);
   
      for j=1:q
         
         a = filter([ coeffMA0(j,i) coeffMA1(j,i) coeffMA2(j,i)],1,u(:,j)) + a;
         
      end
   common(:,i) = a(51:T+50);
end
else
 for i=1:N
   a = zeros(T+50,1);
   
      for j=1:q
         
         b = filter(coeffMA0(j,i),[1 coeffAR1(j,i)],u(:,j));
         b = filter(1,[1 coeffAR2(j,i)],b);
         a = b + a;
         
      end
   common(:,i) = a(51:T+50);   
 end
end
stdcommon = sqrt(sum(var(common))/N);
common = common/stdcommon;
stdidio = sqrt(sum(var(idio))/N);
idio = idio/stdidio;
if opt1 == 'LI' % large idio
idio = idio*1;
elseif opt1 == 'SI' %small idio
    idio = idio*.5;
end
variables = common+idio;



