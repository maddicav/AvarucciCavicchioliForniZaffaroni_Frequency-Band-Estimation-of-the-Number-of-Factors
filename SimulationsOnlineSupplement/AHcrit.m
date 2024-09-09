% kER and KGR determine the number of factors in static approximate factor
% model according to Ahn and Horenstein (2013).
%
%[kER,kGR]=AHcrit(x,qmax)
%
%INPUT:    x          :   T x n data matrix  (required)
%          qmax       :   upper bound on the number of factors (required)  
%         
%OUTPUT:  kER,kGR    :  estimated number of factors as maximizer of ER(k) and GR(k) 
%                            
% -------------------------------------------------------------------------
function[kER,kGR]=AHcrit(x,qmax)
[T,n]=size(x);
% 'Double'-demeaning
M1=mean(x,1);
x=x-M1-mean(x,2)+mean(M1);
%
%
if n<T
S=(x'*x)/(T*n);
else
S=(x*x')/(T*n);
end
%
%
m=min(n,T);
d=eigs(S,m);
%
%
er=d(1:qmax)./d(2:qmax+1);
%
[~,kER]=max(er);
%
%
%
%
sumd=sum(d);
ds1=[sumd sumd-cumsum(d')];
ds2=(ds1(1:qmax))./(ds1(2:qmax+1));
gr=log(ds2(1:(qmax-1)))./(log(ds2(2:qmax)));
%
%
%
[~,kGR]=max(gr);
%
end