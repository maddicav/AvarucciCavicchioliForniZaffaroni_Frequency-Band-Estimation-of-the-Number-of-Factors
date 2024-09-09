%--------------------------------------------------------------------------
%
% Output: TABLE S.6 in Avarucci et al 2024, Online Supplement 
%         Left Panel  : Trend-Cycle Model,  
%         Right Panel : Stop-Band Model
%
%--------------------------------------------------------------------------
%
clear all;
clc;
%
tic
qmax = 8;
n = 120;
T = 240;
%
nrepli = 500;
s=[.6  1.2]; % small and large idiosyncratic shock
kDDR=zeros(6,nrepli,length(s));
ACFZ=zeros(6,3,length(s));
%
for k = 1:length(s)   
for j=1:nrepli
%--------------------------------------------------------------------------
% Stop-Band Model
%--------------------------------------------------------------------------
X =  StopBandModel(n,T,1,s(k));
[kDDR(1,j,k)] =  DDR(X, qmax,.75,[0  0]);       
[kDDR(2,j,k)] =  DDR(X, qmax,.75,[0  2*pi/80]);
[kDDR(3,j,k)] = DDR(X, qmax,.75,[2*pi*20/240 2*pi*21/240]);      
[kDDR(4,j,k)] = DDR(X, qmax,.75,[2*pi/32 2*pi/8]);
[kDDR(5,j,k)] = DDR(X, qmax,.75,[2*pi/8 pi]);
[kDDR(6,j,k)] = DDR(X, qmax,.75,[0 pi]);
% 
[k j]
end   
 ACFZ(:,:,k) = [sum(kDDR(:,:,k)==1,2) sum(kDDR(:,:,k)==2,2) sum(kDDR(:,:,k)>2,2) ]*100/nrepli;

end
BAND = [];
for k=1:length(s)
BAND = [BAND ACFZ(:,:,k)];
end
for k = 1:length(s)   
for j=1:nrepli
%--------------------------------------------------------------------------
% Trend-Cycle Model
%--------------------------------------------------------------------------  
X =  TrendCycleModel(n,T,1,s(k));
[kDDR(1,j,k)] =  DDR(X, qmax,.75,[0  0]);       
[kDDR(2,j,k)] =  DDR(X, qmax,.75,[0  2*pi/80]);
[kDDR(3,j,k)] = DDR(X, qmax,.75,[2*pi*20/240 2*pi*21/240]);      
[kDDR(4,j,k)] = DDR(X, qmax,.75,[2*pi/32 2*pi/8]);
[kDDR(5,j,k)] = DDR(X, qmax,.75,[2*pi/8 pi]);
[kDDR(6,j,k)] = DDR(X, qmax,.75,[0 pi]);
% 
[k j]
end   
 ACFZ(:,:,k) = [sum(kDDR(:,:,k)==1,2) sum(kDDR(:,:,k)==2,2) sum(kDDR(:,:,k)>2,2) ]*100/nrepli;
end

TREND = [];
for k=1:length(s)
TREND = [TREND ACFZ(:,:,k)];
end
%--------------------------------------------------------------------------
%  T A B L E   S.6
%--------------------------------------------------------------------------
%toc
TC = [TREND(:,1:3);TREND(:,4:6)]; %Trend Cycle Model
SB = [BAND(:,1:3);BAND(:,4:6)];     %Stop Band Model
EXP4 = [TC sum(TC,2)  SB sum(SB,2)];
display('Table S.6')
disp(EXP4)
%save EXP4
