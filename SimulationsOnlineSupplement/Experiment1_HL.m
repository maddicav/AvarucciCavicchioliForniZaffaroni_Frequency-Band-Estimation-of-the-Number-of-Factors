%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output: TABLES S.1 and S.2 in Avarucci et al 2024, Online Supplement 
%         Top Panel (MA loadings): set " model='MA' " in line 42
%         Top Panel (AR loadings): set " model='AR' " in line 41
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
%
tic
%--------  SETUP  ---------------------------------------------------------
qmax = 8; 
n = [60 100 70 120 150];
T = [100 100 120 120 120];
m = 15; % Badnwidth for Onatski's test
ln=length(n);
%
nrepli = 500;
%
nq=2; %Different Number of factors (q=2 and q=3) considered in the DGP: [2,3]--> 2-dimensional vector
%--------------------------------------------------------------------------
%
%
% ------ PREALLOCATION ----------------------------------------------------
HL=zeros(nrepli,5,nq);
ON = HL;
ACFZ=zeros(nrepli,5,nq,3);
AH=zeros(nrepli,5,nq,2);
ncorrections=zeros(nrepli,1);
%--------------------------------------------------------------------------
%
% 
%------- SIMULATIONS AND ESTIMATION ---------------------------------------
% Functions required:
% HLmodel -> generates X (T x n)
% ACFZcrit-> estimates 'q' using DER,DGR,DDR (Avarucci et al 2024)
% HLcrit  -> estimates 'q' using Hallin & Liska (2007)
% ONcrit  -> estimates 'q' using Onatski (2009)
% AHcrit  -> estimates 'q' using Ahn & Horenstein (2013)
% choose the model, 'AR' or 'MA' (see the lines below)
%-------------------------------------------------------------------------
model='AR';
% model='MA';
%
for j=1:nrepli
   for i = 1:ln
    for q =2:3
X = HLmodel(model,n(i),T(i),q);
[kDER, kDGR, kDDR,ncorrections(j)] = ACFZcrit(X,qmax);
ACFZ(j,i,(q-1),:) = [kDER, kDGR, kDDR];
%
qlog = HLcrit(X, qmax);
HL(j,i,(q-1)) = qlog;
%
ON(j,i,(q-1))=ONcrit(X,qmax,m);
%
[kER,kGR]=AHcrit(X,qmax);
AH(j,i,(q-1),:) = [kER, kGR];
    end
   end;j
end
%--------------------------------------------------------------------------
%
%
%------- T A B L E  S.1 ---------------------------------------------------   
H_L=[sum(HL(:,:,1)==2,1)'; sum(HL(:,:,2)==3,1)']*100/nrepli;
%
O=[sum(ON(:,:,1)==2,1)'; sum(ON(:,:,2)==3,1)']*100/nrepli;
%
DER=[sum(ACFZ(:,:,1,1)==2,1)'; sum(ACFZ(:,:,2,1)==3,1)']*100/nrepli;
DGR=[sum(ACFZ(:,:,1,2)==2,1)'; sum(ACFZ(:,:,2,2)==3,1)']*100/nrepli;
DDR=[sum(ACFZ(:,:,1,3)==2,1)'; sum(ACFZ(:,:,2,3)==3,1)']*100/nrepli;
%
ER=[sum(AH(:,:,1,1)==2,1)'; sum(AH(:,:,2,1)==3,1)']*100/nrepli;
GR=[sum(AH(:,:,1,2)==2,1)'; sum(AH(:,:,2,2)==3,1)']*100/nrepli;
%
n=n';
T=T';
n=vertcat(n,n);
T=vertcat(T,T);
q=[2*ones(ln,1); 3*ones(ln,1)];
%
%
T1_HL=table(q,n, T, H_L, O, DER, DGR, DDR, ER, GR);
if strcmp(model,'MA')==1
display('Table S.1 MA loadings' )
else
    display('Table S.1 AR loadings' )
    end
disp(T1_HL)

%----- create the xls table -----------------------------------------------
%if strcmp(model,'MA')==1
%filename = 'Tab_HL_MA.xlsx';
%writetable(T1_HL,filename,'Sheet',1)
%else
%    filename = 'Tab_HL_AR.xlsx';
%writetable(T1_HL,filename,'Sheet',1)
%end
%--------------------------------------------------------------------------
%
%
%------- T A B L E  S.2 ---------------------------------------------------  
H_L1=[sum(HL(:,:,1)==1,1)'; sum(HL(:,:,2)==1,1)']*100/nrepli;
H_L2=[sum(HL(:,:,1)==2,1)'; sum(HL(:,:,2)==2,1)']*100/nrepli;
H_L3=[sum(HL(:,:,1)==3,1)'; sum(HL(:,:,2)==3,1)']*100/nrepli;
H_L4=[sum(HL(:,:,1)==4,1)'; sum(HL(:,:,2)==4,1)']*100/nrepli;
H_L5=[sum(HL(:,:,1)>4,1)'; sum(HL(:,:,2)>4,1)']*100/nrepli;
%
O1=[sum(ON(:,:,1)==1,1)'; sum(ON(:,:,2)==1,1)']*100/nrepli;
O2=[sum(ON(:,:,1)==2,1)'; sum(ON(:,:,2)==2,1)']*100/nrepli;
O3=[sum(ON(:,:,1)==3,1)'; sum(ON(:,:,2)==3,1)']*100/nrepli;
O4=[sum(ON(:,:,1)==4,1)'; sum(ON(:,:,2)==4,1)']*100/nrepli;
O5=[sum(ON(:,:,1)>4,1)';  sum(ON(:,:,2)>4,1)']*100/nrepli;
%
DER1=[sum(ACFZ(:,:,1,1)==1,1)'; sum(ACFZ(:,:,2,1)==1,1)']*100/nrepli;
DER2=[sum(ACFZ(:,:,1,1)==2,1)'; sum(ACFZ(:,:,2,1)==2,1)']*100/nrepli;
DER3=[sum(ACFZ(:,:,1,1)==3,1)'; sum(ACFZ(:,:,2,1)==3,1)']*100/nrepli;
DER4=[sum(ACFZ(:,:,1,1)==4,1)'; sum(ACFZ(:,:,2,1)==4,1)']*100/nrepli;
DER5=[sum(ACFZ(:,:,1,1)>4,1)'; sum(ACFZ(:,:,2,1)>4,1)']*100/nrepli;
%
DGR1=[sum(ACFZ(:,:,1,2)==1,1)'; sum(ACFZ(:,:,2,2)==1,1)']*100/nrepli;
DGR2=[sum(ACFZ(:,:,1,2)==2,1)'; sum(ACFZ(:,:,2,2)==2,1)']*100/nrepli;
DGR3=[sum(ACFZ(:,:,1,2)==3,1)'; sum(ACFZ(:,:,2,2)==3,1)']*100/nrepli;
DGR4=[sum(ACFZ(:,:,1,2)==4,1)'; sum(ACFZ(:,:,2,2)==4,1)']*100/nrepli;
DGR5=[sum(ACFZ(:,:,1,2)>4,1)'; sum(ACFZ(:,:,2,2)>4,1)']*100/nrepli;
%
DDR1=[sum(ACFZ(:,:,1,3)==1,1)'; sum(ACFZ(:,:,2,3)==1,1)']*100/nrepli;
DDR2=[sum(ACFZ(:,:,1,3)==2,1)'; sum(ACFZ(:,:,2,3)==2,1)']*100/nrepli;
DDR3=[sum(ACFZ(:,:,1,3)==3,1)'; sum(ACFZ(:,:,2,3)==3,1)']*100/nrepli;
DDR4=[sum(ACFZ(:,:,1,3)==4,1)'; sum(ACFZ(:,:,2,3)==4,1)']*100/nrepli;
DDR5=[sum(ACFZ(:,:,1,3)>4,1)'; sum(ACFZ(:,:,2,3)>4,1)']*100/nrepli;
%
ER1=[sum(AH(:,:,1,1)==1,1)'; sum(AH(:,:,2,1)==1,1)']*100/nrepli;
ER2=[sum(AH(:,:,1,1)==2,1)'; sum(AH(:,:,2,1)==2,1)']*100/nrepli;
ER3=[sum(AH(:,:,1,1)==3,1)'; sum(AH(:,:,2,1)==3,1)']*100/nrepli;
ER4=[sum(AH(:,:,1,1)==4,1)'; sum(AH(:,:,2,1)==4,1)']*100/nrepli;
ER5=[sum(AH(:,:,1,1)>4,1)';  sum(AH(:,:,2,1)>4,1)']*100/nrepli;
%
GR1=[sum(AH(:,:,1,2)==1,1)'; sum(AH(:,:,2,2)==1,1)']*100/nrepli;
GR2=[sum(AH(:,:,1,2)==2,1)'; sum(AH(:,:,2,2)==2,1)']*100/nrepli;
GR3=[sum(AH(:,:,1,2)==3,1)'; sum(AH(:,:,2,2)==3,1)']*100/nrepli;
GR4=[sum(AH(:,:,1,2)==4,1)'; sum(AH(:,:,2,2)==4,1)']*100/nrepli;
GR5=[sum(AH(:,:,1,2)>4,1)';  sum(AH(:,:,2,2)>4,1)']*100/nrepli;
%
T2_HL=table(q,n, T, H_L1, H_L2, H_L3, H_L4, H_L5, O1, O2, O3, O4, O5, ...
    DER1, DER2, DER3, DER4, DER5, DGR1,  DGR2, DGR3, DGR4, DGR5, ...
    DDR1, DDR2, DDR3,DDR4, DDR5, ER1, ER2, ER3, ER4, ER5, GR1, GR2, GR3, GR4, GR5);
if strcmp(model,'MA')==1
display('Table S.2 MA loadings' )
else
    display('Table S.2 AR loadings' )
    end
disp(T2_HL)
%----- create the xls table -----------------------------------------------
%if strcmp(model,'MA')==1
%filename = 'Tab_HL_MA.xlsx';
%writetable(T2_HL,filename,'Sheet',2)
%else
%    filename = 'Tab_HL_AR.xlsx';
%writetable(T2_HL,filename,'Sheet',2)
%end

%toc

