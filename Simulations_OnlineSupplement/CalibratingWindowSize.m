%--------------------------------------------------------------------------
% OUPUT: TABLE S.9 in Avarucci et al 2025, Online Supplement
%
% Requires: DGP3model.m, ACFZcrit.m
%           
%--------------------------------------------------------------------------
clear all;
clc;

tic
qmax = 8;
a = [ .5 .75 1 1.25]; %M=[a*sqrt(T)]
n =  216;
T = [120 160  240];
nrepli = 500;
%
 ACFZ=zeros(nrepli,4,3,3);
%
 for k = 1:3
 for j=1:nrepli
    for i = 1:length(a)
    for q = 1:4
    % X =  modelloARMA(n,T(k),q,'LI','ARMA');
      X=   DGP3model(n,T(k),q,'LI');
    X=standardize(X);
%        
[kER, kGR, kO] = ACFZcrit(X, qmax,a(i));
kDDR = DDR(X, qmax,a(i),[2*pi/32 2*pi/6]);

ACFZ(j,i,q,:) = [kER kGR  kO];
     q
    end
    i
    end
 j
 end

TabSim4bis1(:,:,k)=[   sum(ACFZ(:,:,1,1)==1,1)' sum(ACFZ(:,:,1,2)==1,1)' sum(ACFZ(:,:,1,3)==1,1)';
  sum(ACFZ(:,:,2,1)==2,1)' sum(ACFZ(:,:,2,2)==2,1)' sum(ACFZ(:,:,2,3)==2,1)';
  sum(ACFZ(:,:,3,1)==3,1)' sum(ACFZ(:,:,3,2)==3,1)' sum(ACFZ(:,:,3,3)==3,1)';
sum(ACFZ(:,:,4,1)==4,1)' sum(ACFZ(:,:,4,2)==4,1)' sum(ACFZ(:,:,4,3)==4,1)']*100/nrepli;

 end
display('Table S.9, T=120 (first panel), T=160 (second panel) and T=240 (third panel)')
disp(TabSim4bis1)
%toc


