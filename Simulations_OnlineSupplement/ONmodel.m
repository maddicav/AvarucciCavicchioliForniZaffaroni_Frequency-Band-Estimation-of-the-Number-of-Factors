%--------------------------------------------------------------------------
% [common,idi] = ONmodel(model,N,T,q) simulates the (T x n) common and 
% idiosyncratic components using the DGP described in the second experiment
%
% INPUT    model  :  ('AR' or 'MA')
%          N,T    :  panel dimension
%          q      :  number of factors/shocks
% OUTPUT   idi    :  standardized (T x n) idiosincratic component (see second experiment)
%          common :  standardized (T x n) common component (see second experiment)
%--------------------------------------------------------------------------
function [common,idi] = ONmodel(model,N,T,q)
%--------------------------------------------------------------------------
common = ONcommon(N,T,q,model);
idi = ONidi(N,T);
%--------------------------------------------------------------------------
function idi = ONidi(N,T)
%--------------------------------------------------------------------------
randn('state',sum(100*clock));
rand('state',sum(100*clock));
B=200; %burning sample
rho=rand(N,1)*1.6-0.8;
mdl=arima('Constant',0,'AR',0.2,'Variance',1);
V=simulate(mdl,N,'NumPaths',(T+B));
E=zeros(N,T+B);
for s=2:(T+B)
    E(:,s)=rho.*E(:,(s-1))+V(:,s);
end
E=E(:,B+1:end);
E=E';
idi = (E- ones(T,1)*mean(E))./(ones(T,1)*std(E));
idi=idi(:,randperm(N));
end
%--------------------------------------------------------------------------
function common = ONcommon(N,T,q,model)
%--------------------------------------------------------------------------
randn('state',sum(100*clock));
rand('state',sum(100*clock));
common=zeros(T,N);
    if strcmp(model,'MA')==1
        i2=2;
        factors=randn(T+i2,q);
        coeffMA0=randn(q,N);
        coeffMA1=rand(q,N);
        coeffMA2=rand(q,N);
        for i=1:N
            a = zeros(T+i2,1);
            for j=1:q
                a = coeffMA0(j,i)*filter(conv([1 coeffMA1(j,i)],[1 coeffMA2(j,i)]),1,factors(:,j)) + a;
            end
            common(:,i) = a(i2+1:end);
        end
    end
   %
    if strcmp(model,'AR')==1
        B=200;
        factors=randn(T+B,q);
        coeffAR0=randn(q,N);
        coeffAR1=rand(q,N)*0.1+.8;
        coeffAR2=rand(q,N)*0.1+.5;
        for i=1:N
            a = zeros(T+B,1);
            for j=1:q
                % modified the following line
                a = coeffAR0(j,i)*filter(1,conv([1 -coeffAR1(j,i)],[1 -coeffAR2(j,i)]),factors(:,j)) + a;
            end
            common(:,i) = a(B+1:T+B);
        end
    end
    
    common = (common- ones(T,1)*mean(common))./(ones(T,1)*std(common));
end
end

