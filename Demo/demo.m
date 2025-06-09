clear all
clc
%--------------------------------------------------------------------------
% Part A: Generate the observable data X (T x n)
%--------------------------------------------------------------------------
T=300;       % length of each time series
N=300;        % number of variables 
q=2;         % number of shocks
%--- Generate the Latent Common Component ---------------------------------
burn=100;    % burn-in period 
u = randn(T+burn,q);
b0=randn(q,N);
b1=rand(q,N)*1.6-.8;
b2=rand(q,N)*1.6-.8;
common=zeros(T,N);
%
for i=1:N
   a = zeros(T+burn,1);  
      for j=1:q
         ARfilt = conv([1 -b1(j,i) ],[1 -b2(j,i)]);
         a = b0(j,i)*filter(1,ARfilt,u(:,j)) + a;    
      end
   common(:,i) = a(burn+1:end);
end
%--- Standardization of the Common Component ------------------------------
Normaliz1=((1+b1.*b2).*b0.^2)./((1-b2.^2).*(1-b1.^2).*(1-b1.*b2));% Variance of an AR(2) process
Normaliz2=sum(Normaliz1,1);
scommon=common*diag(Normaliz2.^(-.5));
%--- Generate the Data ----------------------------------------------------
X=scommon+randn(T,N); %The idionsyncratic compoonent is Gaussian iid
%
%
%--------------------------------------------------------------------------
% Part B: Compute DER,DGR and DDR using all Frequencies 
%--------------------------------------------------------------------------
qmax=8; % maximum possible number of shocks 
% qmax must be smaller than min(N,2*M+1)-1, 'M' is the bandwidth
c=.75; % positive constant, M=floor(c*T^(1/2)) 
%
%-- Estimate 'q' by maximizing the DER,DGR and DDR criteria------------------
[kDER, kDGR, kDDR] = ACFZcrit(X, qmax,c); % 'c' optional, default value 0.75
%--------------------------------------------------------------------------
% Part C: Compute DDR on frequency bands 
%--------------------------------------------------------------------------
% Estimate 'q'by maximizing the DDR criteria and using two different frequency bands
freqband1=[0  2*pi/80]; % long-run frequencies
freqband2=[2*pi/32 2*pi/6]; % business-cycle frequencies
kDDR_band1=DDR(X,qmax,c,freqband1); % long-run frequencies
kDDR_band2=DDR(X,qmax,c,freqband2); % business-cycle frequencies
% 'c' is optional, default value 0.75
%--------------------------------------------------------------------------
% Part D: Display the output
%--------------------------------------------------------------------------
disp(['The DER-estimate of the number of shocks using all frequencies is q=', num2str(kDER)]);
disp(' ');
disp(['The DGR-estimate of the number of shocks using all frequencies is q=', num2str(kDGR)]);
disp(' ');
disp(['The DDR-estimate of the number of shocks using all frequencies is q=', num2str(kDDR)]);
disp(' ');
disp(['The DDR-estimate of the number of shocks using long-run frequencies is q=', num2str(kDDR_band1)]);
disp(' ');
disp(['The DDR-estimate of the number of shocks using business-cycle frequencies is q=', num2str(kDDR_band2)]);


