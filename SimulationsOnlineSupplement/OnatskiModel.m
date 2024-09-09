function data = OnatskiModel(N,T,q,sigma2,opt)
% Burn-in period
B=200;
% idiosyncratic terms
%
v = zeros(T+B,N);
e = v;
u = randn(T+B,N);
v(:,1) = u(:,1);
for i = 2:N
    v(:,i) = 0.2*v(:,i-1) + u(:,i);
end
rho = rand(1,N)*1.6-0.8;
e(1,:) = v(1,:);
for t = 2:T+B
    e(t,:) = rho.*e(t-1,:) + v(t,:);
end
e = e((B+1):end,:);
%
% common components
%
F = randn(T+B, q);
if nargin == 4
    opt = 'MA';
elseif nargin == 3
    opt = 'MA';
    sigma2 = 1;
end
%
% model MA
%
common=zeros(T,N);
if opt == 'MA'
    coeffMA0 = randn(q,N);
    coeffMA1 = rand(q,N);
    coeffMA2 = rand(q,N);
for i=1:N,
   a = zeros(T+B,1);
      for j=1:q,
         b = coeffMA0(j,i)*filter(conv([1 coeffMA1(j,i)],[1 coeffMA2(j,i)]),1,F(:,j));
         a = b + a;
      end
   common(:,i) = a((B+1):T+B);   
end
%
% variance normalization
%
common = common*diag(std(common).^(-1)*sqrt(.4+.05*q));
e = sqrt(sigma2)*e*diag(std(e).^(-1)*sqrt(1-(.4+.05*q)));
data = common + e;
end
%
% model AR
%
if opt == 'AR'
    coeffAR0 = randn(q,N);
    coeffAR1 = rand(q,N)*.1 + .8;
    coeffAR2 = rand(q,N)*.1 + .5;
for i=1:N,
   a = zeros(T+B,1);
      for j=1:q,
         b = filter(coeffAR0(j,i),conv([1 coeffAR1(j,i)],[1 coeffAR2(j,i)]),F(:,j));
         a = b + a;
      end
   common(:,i) = a((B+1):T+B);   
end
 %
% variance normalization
%
common = common*diag(std(common).^(-1)*sqrt(.4+.05*q));
e = sqrt(sigma2)*e*diag(std(e).^(-1)*sqrt(1-(.4+.05*q)));
data = common + e;
end
end