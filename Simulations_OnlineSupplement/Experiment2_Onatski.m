%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output: TABLES S.3 and S.4 in Avarucci et al 2024, Online Supplement 
%         Top Panel (MA loadings): set " model='MA' " in line 36
%         Top Panel (AR loadings): set " model='AR' " in line 37
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
%------------------------------------------------------------------------
% General Settings
%------------------------------------------------------------------------
tic
qmax = 8;
n = [70   100  150];
T = [70 120  500];
s2 = [1 2 4 1 2 6 1 8 16];
m = [30  40  65];
h=length(n);
hh=length(s2);
nrepli = 500;
q=2;
%------------------------------------------------------------------------
% Preallocation
%------------------------------------------------------------------------
HL=zeros(nrepli,hh);
ON = HL;
ACFZ=zeros(nrepli,hh,3);
AH=zeros(nrepli,hh,2);
nocorrection=zeros(nrepli,1);
%------- SIMULATIONS AND ESTIMATION ---------------------------------------
% Functions required:
% ONmodel -> generates common (T x n) and idio (T x n), see line 47
% ACFZcrit-> estimates 'q' using DER,DGR,DDR (Avarucci et al 2024)
% HLcrit  -> estimates 'q' using Hallin & Liska (2007)
% ONcrit  -> estimates 'q' using Onatski (2009)
% AHcrit  -> estimates 'q' using Ahn & Horenstein (2013)
% choose the model, 'AR' or 'MA', see the line below
%--------------------------------------------------------------------------
model='MA'; 
%model='AR';
for j=1:nrepli
   for i = 1:h
        [common,idio] =  ONmodel(model,n(i),T(i),q); %generates the common and idiosyncratic components, both with variance 1
        ind=3*(i-1)+1;
      for k=ind:(ind+2)
    X=common+sqrt(s2(k))*idio;
[kER, kGR, kDDR,nocorrection(j)] = ACFZcrit(X, qmax);
ACFZ(j,k,:) = [kER kGR  kDDR];
qlog = HLcrit(X, qmax);
HL(j,k) = qlog;
ON(j,k)=ONcrit(X,qmax,m(i));
%
[kER,kGR]=AHcrit(X,qmax);
AH(j,k,:) = [kER, kGR];
     end
  end;j
end
%------------------------------------------------------------------------
% TABLE S.3
% Top (bottom) panel if model='MA' (model='AR')
%------------------------------------------------------------------------
H_L=(sum(HL(:,:)==q,1)')*100/nrepli;
%
O=(sum(ON(:,:)==q,1)')*100/nrepli;
%
DER=(sum(ACFZ(:,:,1)==q,1)')*100/nrepli;
DGR=(sum(ACFZ(:,:,2)==q,1)')*100/nrepli;
DDR=(sum(ACFZ(:,:,3)==q,1)')*100/nrepli;
%
ER=(sum(AH(:,:,1)==q,1)')*100/nrepli;
GR=(sum(AH(:,:,2)==q,1)')*100/nrepli;
%
nn = [70 70 70  100 100 100 150 150 150]';
TT = [70 70 70 120 120 120 500 500 500]';
%
var=s2';
T1_O=table(nn, TT, var, H_L, O, DER, DGR, DDR, ER, GR);
if strcmp(model,'MA')==1
display('Table S.3 MA loadings' )
else
    display('Table S.3 AR loadings' )
    end
disp(T1_O)
%if strcmp(model,'MA')==1
%filename = 'Tab_Onatski_MA.xlsx';
%writetable(T1_O,filename,'Sheet',1)
%else
%    filename = 'Tab_Onatski_AR.xlsx';
%writetable(T1_O,filename,'Sheet',1)
%end
%------------------------------------------------------------------------
% TABLE S.3
% Top (bottom) panel if model='MA' (model='AR')
%------------------------------------------------------------------------
 H_L1=(sum(HL(:,:)==1,1)')*100/nrepli;
 H_L2=(sum(HL(:,:)==2,1)')*100/nrepli;
 H_L3=(sum(HL(:,:)==3,1)')*100/nrepli;
 H_L4=(sum(HL(:,:)==4,1)')*100/nrepli;
 H_L5=(sum(HL(:,:)>4,1)')*100/nrepli;
%
O1=(sum(ON(:,:)==1,1)')*100/nrepli;
O2=(sum(ON(:,:)==2,1)')*100/nrepli;
O3=(sum(ON(:,:)==4,1)')*100/nrepli;
O4=(sum(ON(:,:)==4,1)')*100/nrepli;
O5=(sum(ON(:,:)>4,1)')*100/nrepli;
%
DER1=(sum(ACFZ(:,:,1)==1,1)')*100/nrepli;
DER2=(sum(ACFZ(:,:,1)==2,1)')*100/nrepli;
DER3=(sum(ACFZ(:,:,1)==3,1)')*100/nrepli;
DER4=(sum(ACFZ(:,:,1)==4,1)')*100/nrepli;
DER5=(sum(ACFZ(:,:,1)>4,1)')*100/nrepli;
%
DGR1=(sum(ACFZ(:,:,2)==1,1)')*100/nrepli;
DGR2=(sum(ACFZ(:,:,2)==2,1)')*100/nrepli;
DGR3=(sum(ACFZ(:,:,2)==3,1)')*100/nrepli;
DGR4=(sum(ACFZ(:,:,2)==4,1)')*100/nrepli;
DGR5=(sum(ACFZ(:,:,2)>4,1)')*100/nrepli;
%% 
DDR1=(sum(ACFZ(:,:,3)==1,1)')*100/nrepli;
DDR2=(sum(ACFZ(:,:,3)==2,1)')*100/nrepli;
DDR3=(sum(ACFZ(:,:,3)==3,1)')*100/nrepli;
DDR4=(sum(ACFZ(:,:,3)==4,1)')*100/nrepli;
DDR5=(sum(ACFZ(:,:,3)>4,1)')*100/nrepli;
%
ER1=(sum(AH(:,:,1)==1,1)')*100/nrepli;
ER2=(sum(AH(:,:,1)==2,1)')*100/nrepli;
ER3=(sum(AH(:,:,1)==3,1)')*100/nrepli;
ER4=(sum(AH(:,:,1)==4,1)')*100/nrepli;
ER5=(sum(AH(:,:,1)>4,1)')*100/nrepli;
%
GR1=(sum(AH(:,:,2)==1,1)')*100/nrepli;
GR2=(sum(AH(:,:,2)==2,1)')*100/nrepli;
GR3=(sum(AH(:,:,2)==3,1)')*100/nrepli;
GR4=(sum(AH(:,:,2)==4,1)')*100/nrepli;
GR5=(sum(AH(:,:,2)>4,1)')*100/nrepli;
%
T2_O=table(nn, TT, var, H_L1, H_L2, H_L3, H_L4, H_L5, O1, O2, O3, O4, O5, ...
    DER1, DER2, DER3, DER4, DER5, DGR1,  DGR2, DGR3, DGR4, DGR5, ...
    DDR1, DDR2, DDR3,DDR4, DDR5, ER1, ER2, ER3, ER4, ER5, GR1, GR2, GR3, GR4, GR5);
if strcmp(model,'MA')==1
display('Table S.4 MA loadings' )
else
    display('Table S.4 AR loadings' )
    end
disp(T2_O)
%
%if strcmp(model,'MA')==1
%filename = 'Tab_Onatski_MA.xlsx';
%writetable(T2_O,filename,'Sheet',2)
%else
%    filename = 'Tab_Onatski_AR.xlsx';
%writetable(T2_O,filename,'Sheet',2)
%end
