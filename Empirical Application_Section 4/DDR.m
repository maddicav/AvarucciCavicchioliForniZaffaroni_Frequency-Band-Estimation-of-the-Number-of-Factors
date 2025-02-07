%kDDR determine the number of factors on a frequency band in GDFMs according to
% Avarucci, Cavicchioli, Forni e Zaffaroni (2025)
%
%[kDDR,DDR,mu11,mum,ncorrections] = DDR(x, qmax,c,band)
%
%INPUT    x          : T x n data matrix  (required)
%         qmax       : upper bound on the number of factors (required)  
%         c          : the bandwidth is computed as M=[c(sqrt(T)) ] (default c=.75)
%         band       : Interval of frequencies considered to compute the estimator
%
%OUTPUT  kDDR    :  estimated number of factors as maximizer of  DDR(k)
%        DDR     :  value of the criteria
%        mu11    :  eingenvalues of the smoothed periodogram computed at
%                   the frequencies [0,M]*2*pi/T
%                   
%        mum          :  eingenvalues of the smoothed periodogram computed at
%                        the frequencies [0,M]*2*pi/T
%        ncorrections :  number of time that the difference of subsequent eigenvalues is smaller than the smallest eigenvalue (denominator DDR)     
%                            
% -------------------------------------------------------------------------


function [kDDR,DDR,mu1,ncorr] = DDR(D,T,qmax,band)
S = floor((T-1)/2);
if nargin < 2
    band = [0 pi];
end

frequencies = (0:2*pi/T:2*pi*S/T)';
bandd = find(frequencies>=band(1) & frequencies<=band(2));

if isempty(bandd)
    disp('the band is empty')
end


for h = 0:S
mu1(h+1,:) = abs(diag(D(:,:,h+1)));
end
if band(1) == 0
weights = [0.5 ones(1,length(bandd)-1)];
weights = weights/sum(weights);
mu = weights*mu1(bandd,:);
mum = weights*mu1(bandd,end);
else
mu = mean(mu1(bandd,:),1);
mum = mean(mu1(bandd,end),1);
end

% compute  DDR and kDDR
den = max((mu(2:qmax+1)-mu(3:qmax+2)),mum);
ncorr = sum(den == mum);
DDR = (mu(1:qmax)-mu(2:qmax+1))./den;
[~,kDDR] = max(DDR);


