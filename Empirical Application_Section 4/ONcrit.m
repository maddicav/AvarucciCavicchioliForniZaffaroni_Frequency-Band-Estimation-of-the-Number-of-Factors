% nfactors estimates the number of factors according to Onatski (2009), Section
% 5.3
%
%nfactors = ONcrit(data,kmax,m)
%
%INPUT:    data       :   T x n data matrix  (required)
%          kmax       :   upper bound on the number of factors (required)  
%          m          :   bandwidth (required)
%
%OUTPUT:  nfactors    :  estimated number of factors 
%         
%                 
%-------------------------------------------------------------------------
function nfactors = ONcrit(data,kmax,m)
data = data';
[n,T] = size(data);
%
MT=floor(.7*sqrt(T));
out=zeros(MT,1);
D=zeros(MT,kmax);
%
%
if nargin==2
    theta=linspace(0,pi-2*pi*(floor(3.3*sqrt(T)+2))/T,MT);
for u=1:MT;
    cval=max(dynamico(data,0,kmax,theta(u)),.01);
for j=1:(kmax-1)
[pval,evalues]=dynamico(data,j,kmax,theta(u));
%
%
%
%
%
        if pval>cval&out(u,1)==0;
            out(u,1)=j;
        else
        end
   D(u,:)=(evalues(1:kmax))'; 
end
if out(u)==0;
    out(u)=kmax;
end
end
%
%
else
theta=linspace(0,pi-2*pi*m/T-.01,MT);
%
%
approx=(1:m)';

for u=1:MT;
    cval=max(dynamico(data,0,kmax,theta(u),approx),.01);
for j=1:kmax-1
[pval,evalues]=dynamico(data,j,kmax,theta(u),approx);
%
%
        if pval>cval&out(u,1)==0;
            out(u,1)=j;
        else
        end
   D(u,:)=(evalues(1:kmax))'; 
end
if out(u)==0;
    out(u)=kmax;
end
end
end
a=(sum(D,2)).^2;
b=min(a);
weights=round(a/b);
C=[out weights];
%[dc,~]=size(C);
c=[];
%
%
for l=1:MT
    cc=C(l,1)*ones(1,C(l,2));
    c=[c,cc];
    nfactors=mode(c);
end
