clc
close all
clear
%-------------------------------------------------------------------------
% Load and transform the variables in FRED-QD (See Section S.5. in Avarucci
% et al 2022)
%-------------------------------------------------------------------------
LoadQ
%
xs = standardize(x(1:end-2,:));
T = size(xs,1);
%--------------------------------------------------------------------------
% Frequency bands of interest
%--------------------------------------------------------------------------
bnd(:,8)=[2*pi/16 2*pi/6];    % short cycles
bnd(:,7)=[2*pi/32 2*pi/16];   % medium cycles
bnd(:,5)=[0 2*pi/80];         % long run
bnd(:,6)=[2*pi/80 2*pi/32];   %long cycles
bnd(:,9)=[2*pi/6 2*pi/2];     % short run
bnd(:,2)=[0 2*pi/6];          %excluding fluctuations less than 18 months
bnd(:,1)=[0 pi];              % all frequencies
bnd(:,4)=[0 0];               % freqency zero
bnd(:,3)=[2*pi/32 2*pi/6];    % business cycle
bndlr=bnd(:,5)';              % long run (lr) frequency band 
bndbc=bnd(:,3)';              % business cycle (bc) frequency band 
qmax=8;
%--------------------------------------------------------------------------
% Estimating the number of shocks using different criteria
%--------------------------------------------------------------------------
[kDER, kDGR, kDDR,~] = ACFZcrit(xs,qmax);
%ACFZ = [kDER, kDGR, kDDR];
qlog = HLcrit(xs, qmax);
HL = qlog;
ON=ONcrit(xs,qmax,20);
%--------------------------------------------------------------------------
% Compute the spectral density matrix and decompose it according to the
% principal component series decomposition (Brillinger 1961)
%--------------------------------------------------------------------------
M = round(.75*sqrt(T));
%M = round(1*sqrt(T));
W = 2*M + 1;
[D,V, Sigma]  = DynamicEigenvalues(xs,W,M);
nof = size(Sigma,3);
for h=1:nof
    for i = 1:size(D,1)
 A(:,:,h,i) = real(V(:,i,h)*D(i,i,h)*V(:,i,h)'); 
    end
end
%--------------------------------------------------------------------------
% Compute DDR on different frequency bands
%--------------------------------------------------------------------------
DD  = DynamicEigenvalues(xs,2*M+1,M);
for j = 1:9
[kDDR(j),DDRALL(:,j),~,ncorr(j)] = DDR(DD,T,qmax,bnd(:,j));
end


kDDR;
S = floor((T-1)/2);
for h = 0:S
mu1(h+1,:) = abs(diag(D(:,:,h+1)));
mum(h+1)=mu1(h+1,end);
end

den = max((mu1(:,2:qmax+1)-mu1(:,3:qmax+2)),mum');
%ncorr = sum(den == mum);
DDRPW = (mu1(:,1:qmax)-mu1(:,2:qmax+1))./den;
[~,kDDRPW] = max(DDRPW,[],2);
cut = floor(nof/2);
%--------------------------------------------------------------------------
% Output: Fig 2
%--------------------------------------------------------------------------
figure(2)
plot(0:2*pi/T:2*pi*cut/T,kDDRPW(1:cut+1,:), 'linewidth',1.5);hold on
plot([ 0.0393 .1374 .6218 2*pi*cut/T],squeeze(kDDR([5  6 3 9])),':o', 'linewidth',2)
limi= axis;
     patch([fliplr(bnd(:,5)') bnd(:,5)'],[limi(3) limi(3) limi(4) limi(4)],[.3 .3 .8]);
    patch([fliplr(bnd(:,3)') bnd(:,3)'],[limi(3) limi(3) limi(4) limi(4)],[.8 .3 .3]);
    %patch([fliplr(bnd(:,8)') bnd(:,8)'],[limi(3) limi(3) limi(4) limi(4)],[.3 .8 .3]);
    alpha 0.1

% rolling sample
k = [1 2 3 ];
passo = 20;
samplesize=160; 
M = round(1*sqrt(samplesize));
for h=1:passo:240-samplesize+1
[kDER0((h-1)/passo+1), kDGR0((h-1)/passo+1)] = ACFZcrit(xs(h:h-1+samplesize,:),qmax,1);
HL0((h-1)/passo+1) = HLcrit(xs(h:h-1+samplesize,:), qmax);
 ON0((h-1)/passo+1)=ONcrit(xs(h:h-1+samplesize,:),qmax,20) ;  
DD  = DynamicEigenvalues(xs(h:h-1+samplesize,:),2*M+1,M);
for j = 1:3
kDDR0((h-1)/passo+1,j) = DDR(DD,samplesize,qmax,bnd(:,k(j)));
end
end
passo = 40;
samplesize=120;
M = round(1*sqrt(samplesize));
for h=1:passo:240-samplesize+1
[kDER1((h-1)/passo+1), kDGR1((h-1)/passo+1)] = ACFZcrit(xs(h:h-1+samplesize,:),qmax,1);
HL1((h-1)/passo+1) = HLcrit(xs(h:h-1+samplesize,:), qmax);
 ON1((h-1)/passo+1)=ONcrit(xs(h:h-1+samplesize,:),qmax,15) ;  
DD  = DynamicEigenvalues(xs(h:h-1+samplesize,:),2*M+1,M);
for j = 1:3
kDDR1((h-1)/passo+1,j) = DDR(DD,samplesize,qmax,bnd(:,k(j)));
end
end
%--------------------------------------------------------------------------
% Monthly Data
%--------------------------------------------------------------------------
clear x
LoadM2
mxs = standardize(x(1:end,:));

mT = size(mxs,1);
% the bands of interest


% different criteria
[mkDER, mkDGR,~,~] = ACFZcrit(mxs,qmax);
mqlog = HLcrit(mxs, qmax);
mHL = mqlog;
mON=ONcrit(mxs,qmax,20);

% DDR on different frequency bands 

    
 MM = round(.75*sqrt(mT));
mDD  = DynamicEigenvalues(mxs,2*MM+1,MM);
for j = 1:3
[mkDDR(j),~,~,] = DDR(mDD,mT,qmax,bnd(:,j));
end

%--------------------------------------------------------------------------
% Output: Table 1
%--------------------------------------------------------------------------
Tab1=[kDDR([1 2 3]) kDER kDGR HL ON; kDDR0 kDER0' kDGR0' HL0' ON0'; kDDR1 kDER1' kDGR1' HL1' ON1'; mkDDR([1 2 3]) mkDER mkDGR mHL mON];

rowtitle{1} = '1960Q2-2020Q1';   
rowtitle{2} = '1960Q2-2000Q1';  
rowtitle{3} = '1965Q2-2005Q1';  
rowtitle{4} = '1970Q2-2010Q1';  
rowtitle{5} = '1975Q2-2015Q1';  
rowtitle{6} = '1980Q2-2020Q1'; 
rowtitle{7} = '1960Q2-1990Q1'; 
rowtitle{8} = '1970Q2-2000Q1'; 
rowtitle{9} = '1980Q2-2010Q1'; 
rowtitle{10} = '1990Q2-2020Q1'; 
rowtitle{11} = '1960M1-2020M3'; 
ctit{1} = 'DDR';
ctit{2} = 'DDRa';
ctit{3} = 'DDRbc';
ctit{4} = 'DER';
ctit{5} = 'DGR';
ctit{6} = 'HL';
ctit{7} = 'O';
%--------------------------------------------------------------------------
% Output: Table 1
%--------------------------------------------------------------------------
disp('TABLE 1')
MakeTable(Tab1,0,rowtitle,ctit);
%--------------------------------------------------------------------------
% Output: Figure 3
%--------------------------------------------------------------------------
k = [1 2 3 ];
passo = 1;
samplesize=160; 
M = round(1*sqrt(samplesize));
for h=1:passo:240-samplesize+1
% [kDER0((h-1)/passo+1), kDGR0((h-1)/passo+1)] = CFcriterion(xs(h:h-1+samplesize,:),qmax,1);
% [~, HL0((h-1)/passo+1)] = HallinLiska(xs(h:h-1+samplesize,:), qmax);
%  AA0((h-1)/passo+1)=Alexei2(xs(h:h-1+samplesize,:),qmax,20) ;  
DD  = DynamicEigenvalues(xs(h:h-1+samplesize,:),2*M+1,M);
for j = 1:3
kDDRrol((h-1)/passo+1,j) = DDR(DD,samplesize,qmax,bnd(:,k(j)));
end
end
figure(3) 
plot([1960+.375:.25:1960+.375+(240-samplesize)*.25], kDDRrol(:,1),'go','MarkerSize',12,'LineWidth',2);hold on
plot([1960+.375:.25:1960+.375+(240-samplesize)*.25], kDDRrol(:,2),'bo','MarkerSize',8,'LineWidth',1.5);
plot([1960+.375:.25:1960+.375+(240-samplesize)*.25], kDDRrol(:,3),'r+','LineWidth',1.5);
axis([1960 1981 0.5 3.5]);grid
%--------------------------------------------------------------------------
% Output: Table 2
%--------------------------------------------------------------------------
S = nof-1;
for h=1:nof
sig(:,h) = diag(Sigma(:,:,h));     
mu(:,h) = diag(D(:,:,h));

end
gdp1 = squeeze(A(1,1,:,:))';
cons1 = squeeze(A(2,2,:,:))';
inv1 = squeeze(A(7,7,:,:))';
dunrate1 = squeeze(A(57,57,:,:))';
hours1 = squeeze(A(69,69,:,:))';
gdpdef1 = squeeze(A(86,86,:,:))';
ffr1 = squeeze(A(128,128,:,:))';
frequencies = (0:2*pi/T:2*pi*S/T)';
for i = 1 : 8
band = find(frequencies>=bnd(1,i) & frequencies<=bnd(2,i));
eigen(:,i) = mean(mu(:,band),2);
gdp(:,i) = mean(gdp1(:,band),2);
cons(:,i) = mean(cons1(:,band),2);
inv(:,i) = mean(inv1(:,band),2);
du(:,i) = mean(dunrate1(:,band),2);
hours(:,i) = mean(hours1(:,band),2);
gdpdef(:,i) = mean(gdpdef1(:,band),2);
ffr(:,i) = mean(ffr1(:,band),2);
end
vardect = eigen(1:6,[1  5:8 3])*100./sum(eigen(:,[1  5:8 3]));
vardecgdp = gdp(1:6,[1  5:8 3])*100./sum(gdp(:,[1  5:8 3]));
vardeccons = cons(1:6,[1  5:8 3])*100./sum(cons(:,[1  5:8 3]));
vardecinv = inv(1:6,[1  5:8 3])*100./sum(inv(:,[1  5:8 3]));

vardecdu = du(1:6,[1  5:8 3])*100./sum(du(:,[1  5:8 3]));
vardechours = hours(1:6,[1  5:8 3])*100./sum(hours(:,[1  5:8 3]));
vardecgdpdef = gdpdef(1:6,[1  5:8 3])*100./sum(gdpdef(:,[1  5:8 3]));
vardecffr = ffr(1:6,[1  5:8 3])*100./sum(ffr(:,[1  5:8 3]));
% total explained variance
%Tab3 = vardect
ctit{1} = 'All frequencies';
ctit{2} = 'Long run';
ctit{3} = 'Long cycles';
ctit{4} = 'Medium cycles';
ctit{5} = 'Short cycles';
ctit{6} = 'Business cycle';
%rtit{1} = '1st dynamic PC';
%rtit{2} = '2nd dynamic PC';
%rtit{3} = '3rd dynamic PC';
%rtit{4} = '4th dynamic PC';
%rtit{5} = '5th dynamic PC';
%rtit{6} = '6th dynamic PC';
%MakeTable(Tab3,1,rtit,ctit);
% how much we loose of some key macro variables by ignoring the third, the
% fourth and the fifth dpcs
Tab2 = [sum(vardecgdp(1:2,:));sum(vardecgdp(3:5,:));sum(vardeccons(1:2,:));sum(vardeccons(3:5,:));...
sum(vardecinv(1:2,:));sum(vardecinv(3:5,:));sum(vardecdu(1:2,:));sum(vardecdu(3:5,:));...
sum(vardechours(1:2,:));sum(vardechours(3:5,:));...
    sum(vardecgdpdef(1:2,:));sum(vardecgdpdef(3:5,:)); sum(vardecffr(1:2,:));sum(vardecffr(3:5,:))];
axis([1960 1981 0.5 3.5]);grid
%--------------------------------------------------------------------------
% Output Table 3
%--------------------------------------------------------------------------
rtit{1} = 'first 2';
rtit{2} = 'next 3';
rtit{3} = 'first 2';
rtit{4} = 'next 3';
rtit{5} = 'first 2';
rtit{6} = 'next 3';
rtit{7} = 'first 2';
rtit{8} = 'next 3';
rtit{9} = 'first 2';
rtit{10} = 'next 3';
rtit{11} = 'first 2';
rtit{12} = 'next 3';
rtit{13} = 'first 2';
rtit{14} = 'next 3';
disp('TABLE 2')
MakeTable(Tab2,1,rtit,ctit);

% variance explained by the first and the second principal component for a
% few key macrovariables
Tab3 = [vardecgdp(1:2,:);vardeccons(1:2,:);...
    vardecinv(1:2,:);vardecdu(1:2,:);...
    vardechours(1:2,:);vardecgdpdef(1:2,:);...
    vardecffr(1:2,:)];
rtit{1} = '1st';
rtit{2} = '2nd';
rtit{3} = '1st';
rtit{4} = '2nd';
rtit{5} = '1st';
rtit{6} = '2nd';
rtit{7} = '1st';
rtit{8} = '2nd';
rtit{9} = '1st';
rtit{10} = '2nd';
rtit{11} = '1st';
rtit{12} = '2nd';
rtit{13} = '1st';
rtit{14} = '2nd';
display('TABLE 3')
MakeTable(Tab3,1,rtit,ctit);
axis([1960 1981 0.5 3.5]);grid
%--------------------------------------------------------------------------
% OutpuT: Figure 4
%--------------------------------------------------------------------------
tit{1} ='GDP';
tit{2} ='Consumption';
tit{3} ='Investment';
tit{4} ='Unemployment';
tit{5} ='Hours worked';
tit{6} ='Inflation';
tit{7} ='Federal funds rate';
tit{8} ='GDP-Inflation cospectrum';
cut = floor(nof/2);
vari = [1 2 7 57 69 86 128];

figure(4)
colors = [0 0 0; 1 0 0 ; 0 0 1];
colororder(colors)
for ss = 1:4
    subplot(2,2,ss)
Ass = squeeze(A(vari(ss),vari(ss),:,:));

plot(0:2*pi/T:2*pi*cut/T,[squeeze(Sigma(vari(ss),vari(ss),1:cut+1)) Ass(1:cut+1,1:2) ], ...
    'Linewidth',1)
hold on
axis tight
limi= axis;
     patch([fliplr(bndlr) bndlr],[limi(3) limi(3) limi(4) limi(4)],[.3 .3 .7]);
    patch([fliplr(bndbc) bndbc],[limi(3) limi(3) limi(4) limi(4)],[.7 .3 .3]);
    alpha 0.1
 title(tit{ss})  
grid 'minor'
end




for h=1:nof
    for i=1:2
B(i,h) = real(V(1,i,h)*D(i,i,h)*(1-exp(1i*2*pi*h/T))*V(86,i,h)');

    end
end

figure(5)
colororder(colors)
for ss = 1:3
    subplot(2,2,ss)
Ass = squeeze(A(vari(4+ss),vari(4+ss),:,:));
plot(0:2*pi/T:2*pi*cut/T,[squeeze(Sigma(vari(4+ss),vari(4+ss),1:cut+1)) Ass(1:cut+1,1:2) ],'Linewidth',1)
hold on
axis tight
limi= axis;
     patch([fliplr(bndlr) bndlr],[limi(3) limi(3) limi(4) limi(4)],[.3 .3 .7]);
    patch([fliplr(bndbc) bndbc],[limi(3) limi(3) limi(4) limi(4)],[.7 .3 .3]);
    alpha 0.1
  title(tit{ss+4})   
grid 'minor'
end
subplot(2,2,4)
plot(0:2*pi/T:2*pi*cut/T,B(1,1:cut+1)','Linewidth',1,'color',[1 0 0]); hold on
plot(0:2*pi/T:2*pi*cut/T,B(2,1:cut+1)','y','Linewidth',1,'color',[0 0 1])
axis tight
limi= axis;
     patch([fliplr(bndlr) bndlr],[limi(3) limi(3) limi(4) limi(4)],[.3 .3 .7]);
    patch([fliplr(bndbc) bndbc],[limi(3) limi(3) limi(4) limi(4)],[.7 .3 .3]);
    alpha 0.1
   title(tit{8})  
grid 'minor'

