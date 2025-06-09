%kDDR determine the number of factors on a frequency band in GDFMs according to
% Avarucci, Cavicchioli, Forni e Zaffaroni (2025)
%
%[kDDR,DDR,mu11,mum,ncorrections] = DDR(x, qmax,c,band)
%
%INPUTS   x          : T x N data matrix  (required)
%                      The time series dimension 'T' equals the number of rows.
%                      The cross-sectional dimension 'N' equals the number of columns.
%         qmax       : upper bound on the number of factors (required)  
%         c          : the bandwidth is computed as M=[c(sqrt(T)) ] (default c=.75)
%         band       : Interval of frequencies considered to compute the estimator
%
%OUTPUTS  kDDR    :  estimated number of factors using DDR(k) criterion
%         DDR     :  value of the criteria
%         mu11    :  eingenvalues of the smoothed periodogram computed at
%                   the frequencies [0,M]*2*pi/T
%                   
%         mum          :  eingenvalues of the smoothed periodogram computed at
%                        the frequencies [0,M]*2*pi/T
%         ncorrections :  number of time that the difference of subsequent eigenvalues is smaller than the smallest eigenvalue (denominator DDR)     
%                            
% -------------------------------------------------------------------------

function [kDDR,DDR,mu11,mum,ncorrections] = DDR(x, qmax,c,band)
[T, ~] = size(x);
S = floor((T-1)/2);
if nargin < 4
    band = [0 pi];
end
frequencies = (0:2*pi/T:2*pi*S/T)';
bandd = find(frequencies>=band(1) & frequencies<=band(2));
if isempty(bandd)
    disp('the band is empty')
end
if nargin < 3
c = .75;
end
[T,~] = size(x);
M = round(sqrt(T)*c);
if nargin <2
    qmax = 2*M+1;
end

[mu11,mum]  = DynamicEigenvaluesPERS(x,qmax,M);

mu = mean(mu11(bandd,:),1);
mum = mean(mum(bandd,:),1);

% compute  DDR and kDDR
den = max((mu(2:end-1)-mu(3:end)),mum);
ncorrections = sum(den == mum);
DDR = (mu(1:end-2)-mu(2:end-1))./den;
[~,kDDR] = max(DDR);

