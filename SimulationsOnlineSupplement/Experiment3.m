%-----------------------------------------------------------------------------------
%
% Output: TABLE S.5 in Avarucci et al 2024, Online Supplement 
%         Top Panel (large idiosyncratic compoment): set " sizeidio='LI' " - Line 32
%         Top Panel (small idiosyncratic compoment): set " sizeidio='SI' " - Line 33
%
%------------------------------------------------------------------------------------

clear all;
clc;

tic
qmax = 8;
n = [60 120 60 120 240];
T = [120 80 240 240 480];
m = [15 15 20 20 30]; 
nrepli = 500;
HL=zeros(nrepli,5,3);
ON = HL;
ACFZ=zeros(nrepli,5,3,3);
AH=zeros(nrepli,5,3,2);
%------- SIMULATIONS AND ESTIMATION ---------------------------------------
% Functions required:
% DGP3model -> generates X (T x n)
% ACFZcrit-> estimates 'q' using DER,DGR,DDR (Avarucci et al 2024)
% HLcrit  -> estimates 'q' using Hallin & Liska (2007)
% ONcrit  -> estimates 'q' using Onatski (2009)
% AHcrit  -> estimates 'q' using Ahn & Horenstein (2013)
% choose the relative variance of the idiosyncratic component (sizeidio), 'LI' or 'SI' (see the lines below)
%-------------------------------------------------------------------------
%sizeidio='LI'; % Large idiosyncratic component
sizeidio='SI'; % Small idiosyncratic component
for j=1:nrepli
   for i = 1:length(n)       
    for q = 1:3
     X =  DGP3model(n(i),T(i),q*2,sizeidio);
     X=standardize(X);
%       
[kER, kGR, kDDR] = ACFZcrit(X, qmax);
ACFZ(j,i,q,:) = [kER kGR  kDDR];
%
[kER,kGR]=AHcrit(X,qmax);
AH(j,i,q,:) = [kER, kGR];
%
qlog = HLcrit(X, qmax);
HL(j,i,q) = qlog;
%
ON(j,i,q)=ONcrit(X,8,m(i));
    end
   end;j
end
%
%
%
H_Ls=[sum(HL(:,:,1)<2,1)';sum(HL(:,:,2)<4,1)';sum(HL(:,:,3)<6,1)']*100/nrepli;   % underestimated q
H_Lc=[sum(HL(:,:,1)==2,1)';sum(HL(:,:,2)==4,1)';sum(HL(:,:,3)==6,1)']*100/nrepli; % correctly estimated q
H_Lb=[sum(HL(:,:,1)>2,1)';sum(HL(:,:,2)>4,1)';sum(HL(:,:,3)>6,1)']*100/nrepli;   % overestimated q
%
Os=[sum(ON(:,:,1)<2,1)';sum(ON(:,:,2)<4,1)';sum(ON(:,:,3)<6,1)']*100/nrepli;
Oc=[sum(ON(:,:,1)==2,1)';sum(ON(:,:,2)==4,1)';sum(ON(:,:,3)==6,1)']*100/nrepli;
Ob=[sum(ON(:,:,1)>2,1)';sum(ON(:,:,2)>4,1)';sum(ON(:,:,3)>6,1)']*100/nrepli;
%
DERs=[sum(ACFZ(:,:,1,1)<2,1)';sum(ACFZ(:,:,2,1)<4,1)';sum(ACFZ(:,:,3,1)<6,1)']*100/nrepli;
DERc=[sum(ACFZ(:,:,1,1)==2,1)';sum(ACFZ(:,:,2,1)==4,1)';sum(ACFZ(:,:,3,1)==6,1)']*100/nrepli;
DERb=[sum(ACFZ(:,:,1,1)>2,1)';sum(ACFZ(:,:,2,1)>4,1)';sum(ACFZ(:,:,3,1)>6,1)']*100/nrepli;
%
DGRs=[sum(ACFZ(:,:,1,2)<2,1)';sum(ACFZ(:,:,2,2)<4,1)';sum(ACFZ(:,:,3,2)<6,1)']*100/nrepli;
DGRc=[sum(ACFZ(:,:,1,2)==2,1)';sum(ACFZ(:,:,2,2)==4,1)';sum(ACFZ(:,:,3,2)==6,1)']*100/nrepli;
DGRb=[sum(ACFZ(:,:,1,2)>2,1)';sum(ACFZ(:,:,2,2)>4,1)';sum(ACFZ(:,:,3,2)>6,1)']*100/nrepli;
%
DDRs=[sum(ACFZ(:,:,1,3)<2,1)';sum(ACFZ(:,:,2,3)<4,1)';sum(ACFZ(:,:,3,3)<6,1)']*100/nrepli;
DDRc=[sum(ACFZ(:,:,1,3)==2,1)';sum(ACFZ(:,:,2,3)==4,1)';sum(ACFZ(:,:,3,3)==6,1)']*100/nrepli;
DDRb=[sum(ACFZ(:,:,1,3)>2,1)';sum(ACFZ(:,:,2,3)>4,1)';sum(ACFZ(:,:,3,3)>6,1)']*100/nrepli;
%
ERs=[sum(AH(:,:,1,1)<2,1)';sum(AH(:,:,2,1)<4,1)';sum(AH(:,:,3,1)<6,1)']*100/nrepli;
ERc=[sum(AH(:,:,1,1)==2,1)';sum(AH(:,:,2,1)==4,1)';sum(AH(:,:,3,1)==6,1)']*100/nrepli;
ERb=[sum(AH(:,:,1,1)>2,1)';sum(AH(:,:,2,1)>4,1)';sum(AH(:,:,3,1)>6,1)']*100/nrepli;
%
GRs=[sum(AH(:,:,1,2)<2,1)';sum(AH(:,:,2,2)<4,1)';sum(AH(:,:,3,2)<6,1)']*100/nrepli;
GRc=[sum(AH(:,:,1,2)==2,1)';sum(AH(:,:,2,2)==4,1)';sum(AH(:,:,3,2)==6,1)']*100/nrepli;
GRb=[sum(AH(:,:,1,2)>2,1)';sum(AH(:,:,2,2)>4,1)';sum(AH(:,:,3,2)>6,1)']*100/nrepli;
%
n=n';
T=T';
nn=vertcat(n,n,n);
TT=vertcat(T,T,T);
ln=length(n);
qq=[2*ones(ln,1); 4*ones(ln,1); 6*ones(ln,1) ];
TT_EXP3=table(qq,nn, TT, H_Ls, H_Lc, H_Lb,  Os, Oc, Ob, DERs, DERc, DERb,  DGRs, DGRc,  DGRb, DDRs, DDRc, DDRb,  ERs, ERc, ERb,  GRs, GRc, GRb );
if strcmp(sizeidio,'SI')==1; 
display('Table S.5 Small idiosyncratic component' )
else
    display('Table S.5 Large idiosyncratic component' )
    end
disp(TT_EXP3)
%if strcmp(sizeidio,'SI')==1
%filename = 'Tab_EXP3_small.xlsx';
%writetable(TT_EXP3,filename,'Sheet',2)
%else
%    filename = 'Tab_EXP3_large.xlsx';
%writetable(TT_EXP3,filename,'Sheet',2)
%end
%
%
%toc