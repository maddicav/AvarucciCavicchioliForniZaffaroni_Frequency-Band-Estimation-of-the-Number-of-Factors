%-------------------------------------------------------------------------
% panel = HLmodel(model,N,T,q) simulates a (T x n) panel data using the DGP
% described in the first experiment
%
% INPUT    model  :  ('AR' or 'MA')
%          N,T    :  panel dimension
%          q      :  number of factors/shocks
% OUTPUT   panel  : (T x n) panel data (see first experiment)
%          common : (T x n) common component in panel (see first experiment)
%-------------------------------------------------------------------------
function panel = HLmodel(model,N,T,q)
%-------------------------------------------------------------------------
% simulate the (T x n) idiosyncratic component
idi = xidi(N,T);
idi = (idi- ones(T,1)*mean(idi))./(ones(T,1)*std(idi));
idi = idi(:,randperm(N));
% simulate the (T x n) idiosyncratic component
common = xcom(N,T,q,model);
common = (common- ones(T,1)*mean(common))./(ones(T,1)*std(common));
%
panel = sqrt(0.5)*common + sqrt(0.5)*idi;  
 
%-------------------------------------------------------------------------
function idi = xidi(N,T)
%-------------------------------------------------------------------------
%
i2=2;
e1=randn(T+i2,N+4);
e2=e1(1+i2:T+i2,:);

C1 = diag(ones(N,1),0)+diag(ones(N-1,1),1)+diag(ones(N-2,1),2)+diag(ones(N-3,1),3)+diag(ones(N-4,1),4);
A=ones(4,4);
A=tril(A);
B=zeros(4,N-4);
C1=[A B; C1];
C2 = C1;

for i1=1:i2
    e2=[e2 e1(1+i2-i1:T+i2-i1,:)];
    C2 = [C2; C1];
end

coeff = (0.5 + rand((i2+1)*(N+4),N)).*C2;
idi = e2*coeff;   
idi = (idi- ones(T,1)*mean(idi))./(ones(T,1)*std(idi));
end
%-------------------------------------------------------------------------
function common = xcom(N,T,q,model)
%-------------------------------------------------------------------------
% q must be equal or smaller than 3!
randn('state',sum(100*clock));
rand('state',sum(100*clock));
%
common=zeros(T,N);
    if strcmp(model,'MA')==1
        i2=2;
        %factors=randn(T+i2,q);
        D = [1,0.5, 1.5];
        factors=randn(T+i2,q).*(ones(T+i2,1)*sqrt(D(1:q)));
        sfactors=factors(1+i2:T+i2,:);
        if i2~=0
            for i1=1:i2
                sfactors=[sfactors factors(1+i2-i1:T+i2-i1,:)];
            end
        end
        coeff=randn((i2+1)*q,N);
        common=sfactors*coeff;
%    
    elseif strcmp(model,'AR')==1
        D = [1,0.5, 1.5];
        factors=randn(T+50,q).*(ones(T+50,1)*sqrt(D(1:q)));
        % factors=randn(T+50,q);
        coeff1=randn(q,N);
        coeff2=rand(q,N)*0.1+.8;% 0.75
        coeff3=rand(q,N)*0.1+.5;% 0.75
        for i=1:N
            a = zeros(T+50,1);
            for j=1:q
                % modified the following line
                a = coeff1(j,i)*filter(1,[1 -1*coeff3(j,i)],filter(1,[1 -1*coeff2(j,i)],factors(:,j))) + a;
            end
            common(:,i) = a(51:T+50);
        end
    else error('please choose the AR or MA model')
    end
    end
end


