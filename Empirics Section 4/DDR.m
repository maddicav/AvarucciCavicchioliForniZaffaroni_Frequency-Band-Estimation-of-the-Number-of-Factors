% ESTIMATION WITH PERIODOGRAM SMOOTHING

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


