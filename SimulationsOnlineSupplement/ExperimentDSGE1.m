%-----------------------------------------------------------------------------------
%
% Output: TABLE S.7 in Avarucci et al 2024, Online Supplement 
%        
% The files results_jpt.mat and solution_acd have been kindly provided by
% Fabrice Collard and are available at http://fabcol.free.fr/research.html
% under "Replication Material - Business Cycle Anatomy".
%------------------------------------------------------------------------------------
clear 
close all 
clc
tic
TT =  [120 240 120 240];
nn = [ 60 120 60 120];
M = .75;
qq = [4 4 7 7];

qmax = 15;
nrepli = 500;
ss = [.1 .2];

%% JPT
load results_jpt;

diaJPT=diag(sqrtm(sol.Sig));

for h = 1:2
    s = ss(h);
for j = 1:length(nn)
    n=nn(j);
    T = TT(j);
     q = qq(j);
     dia = diaJPT;
     dia(q+1:end)=0;
sol.Sig = diag(dia);
for rep = 1:nrepli

t=1;

    
DATA=zeros(11,T);
X=zeros(22,T);
eps=randn(7,T);

X(:,1)=sol.Me*sol.Sig*eps(:,1);
DATA(:,1)=sol.My*X(:,1);
for t=2:T
    X(:,t)=sol.Mx*X(:,t-1)+sol.Me*sol.Sig*eps(:,t);
    DATA(:,t)=sol.My*X(:,t);
end
x = DATA';
m = size(x,2);
x(2:end,[1 2 3 5 6 11    ]) = diff((x(:,[1 2 3 5 6 11     ])));
x = standardize(x(2:end,:));
coeff = rand(m+q,n-m)*2-1;
%
xadd = [x eps(1:q,2:end)']*coeff ;
xx = [x xadd];
e = randn(T,n);
errors = e(2:T,:);
for jj=1:n
errors(:,jj)= filter(1,[1 -.5],e(2:T,jj));
end
errors=standardize(errors)*sqrt(s);
xs = standardize(xx)+ errors;
xs = standardize(xs);
[kDER(rep,j,h), kDGR(rep,j,h), kDDR(rep,j,h)] = ACFZcrit(xs, qmax,M);
%qlog = HLcrit(xs, qmax);
%HL(rep,j,h) = qlog;
rep
end
end
end
%HLq1 = sum([HL(:,1:4,1) HL(:,1:4,2)]==qq(1),1)';
kDDRq1 = sum([kDDR(:,1:2,1) kDDR(:,1:2,2)]==qq(1),1)';
kDERq1 = sum([kDER(:,1:2,1) kDER(:,1:2,2)]==qq(1),1)';
kDGRq1 = sum([kDGR(:,1:2,1) kDGR(:,1:2,2)]==qq(1),1)';
r1 = [kDERq1 kDGRq1 kDDRq1];
r1 = [r1(1:2,: ) r1(3:4,: )];
%HLq2 = sum([HL(:,5:8,1) HL(:,5:8,2)]==qq(5),1)';
kDDRq2 = sum([kDDR(:,3:4,1) kDDR(:,3:4,2)]==qq(3),1)';
kDERq2 = sum([kDER(:,3:4,1) kDER(:,3:4,2)]==qq(3),1)';
kDGRq2 = sum([kDGR(:,3:4,1) kDGR(:,3:4,2)]==qq(3),1)';
r2 = [kDERq2 kDGRq2 kDDRq2];
r2 = [r2(1:2,: ) r2(3:4,: )];
resultJPT = [r1;r2]*100/nrepli;

clear HL kDER kDDR kDGR
%% ACD

qq = [4 4 8 8];
load solution_acd;
diaacd=diag(sqrtm(sol.sig));

for h = 1:2
    s = ss(h);
for j = 1:length(nn)
    n=nn(j);
    T=TT(j);
    q=qq(j);
    dia=diaacd;
    dia(q+1:end)=0;
    sol.sig = diag(dia);
for rep = 1:nrepli
DATA=zeros(10,T);
inf=zeros(T,1);
X=zeros(11,T);
eps=randn(8,T);
epsinf=randn(1,T);
inf(1)=0.27*epsinf(1);
X(:,1)=sol.me*(sol.sig)*eps(:,1);
DATA(:,1)=sol.my*X(:,1);
for t=2:T
    X(:,t)=sol.mx*X(:,t-1)+sol.me*(sol.sig)*eps(:,t);
    DATA(:,t)=sol.my*X(:,t);
    inf(t)=0.89*inf(t-1)+0.27*epsinf(t);


end

%DATAwithI=[DATA ; inf'];
%eps = [eps;epsinf];

x = DATA';
m = size(x,2);
x(2:end,[1 2 3  5 6 8 10]) = diff((x(:,[1 2 3  5 6 8 10])));
x = standardize(x(2:end,:));
coeff = rand(m+q,n-m)*2-1;
xadd = [x eps(1:q,2:end)']*coeff ;
xx = [x xadd];
e = randn(T,n);
errors = e(2:T,:);
for jj=1:n
errors(:,jj)= filter(1,[1 -.5],e(2:T,jj));
end
errors=standardize(errors)*sqrt(s);
xs = standardize(xx)+ errors;
xs = standardize(xs);
[kDER(rep,j,h), kDGR(rep,j,h), kDDR(rep,j,h)] = ACFZcrit(xs, qmax,M);
rep
end
end
end
kDDRq1 = sum([kDDR(:,1:2,1) kDDR(:,1:2,2)]==qq(1),1)';
kDERq1 = sum([kDER(:,1:2,1) kDER(:,1:2,2)]==qq(1),1)';
kDGRq1 = sum([kDGR(:,1:2,1) kDGR(:,1:2,2)]==qq(1),1)';
r1 = [kDERq1 kDGRq1 kDDRq1];
r1 = [r1(1:2,: ) r1(3:4,: )];
kDDRq2 = sum([kDDR(:,3:4,1) kDDR(:,3:4,2)]==qq(3),1)';
kDERq2 = sum([kDER(:,3:4,1) kDER(:,3:4,2)]==qq(3),1)';
kDGRq2 = sum([kDGR(:,3:4,1) kDGR(:,3:4,2)]==qq(3),1)';
r2 = [kDERq2 kDGRq2 kDDRq2];
r2 = [r2(1:2,: ) r2(3:4,: )];
resultACD = [r1;r2]*100/nrepli;



result = [ resultJPT ; resultACD];
display('Table S.7')
disp(result)



