% kDER, KDGR, and kDDR determine the number of factors (including all frequencies) in GDFMs according to
% Avarucci, Cavicchioli, Forni e Zaffaroni (2025)
%
%[kDER, kDGR, kDDR,ncorrections,mu] = ACFZcrit(x, qmax,c)
%
%INPUTS    x          :   T x n data matrix  (required)
%          qmax       :   upper bound on the number of factors (required)  
%          c          :   the bandwidth is computeas as M=[c(sqrt(T)) ] (default c=.75)
%
%OUTPUTS  kDER,kDGR,kDDR    :  estimated number of factors as maximizer of DER(k), DGR(k) and DDR(k)
%         ncorrections      :  number of time that the difference of subsequent eigenvalues is smaller than the smallest eigenvalue (denominator DDR) 
%         mu                :  average (across frequencies) of the 2M+1 (see line 25) largest smoothed periodogram eigenvalues. 
%  ------------------------------------------------------------------------                          
% It requires the function "DynamicEigenvaluesPERS.m" to compute the
% eigenvalues of the smoothed periodgram
% -------------------------------------------------------------------------

function [kDER, kDGR, kDDR,ncorrections,mu] = ACFZcrit(x, qmax,c)


if nargin < 3
c = .75;
end
[T,~] = size(x);
M = round(c*sqrt(T));
if nargin <2
    qmax = 2*M+1;
end
%--------------------------------------------------------------------------
% Compute the (2M+1) largest eigenvalues (columns) of the smoothed periodogram at the
% Fourier frequencies (rowd)
%
[mu11,mum]  = DynamicEigenvaluesPERS(x,(2*M+1),M);
%--------------------------------------------------------------------------

% compute DER and kDER
mu = mean(mu11,1);
mum = mean(mum);
DER = mu(1:qmax)./mu(2:(qmax+1));
[~,kDER] = max(DER);
%kER=kER-1;

% compute DGR and kDGR
%V = V+mu(1);
V = sum(mean(mu11));
dd = [V V - cumsum(mu)];
mustar = dd(1:(qmax))./dd(2:(qmax+1));
DGR = log(mustar(1:(qmax-1)))./log(mustar(2:(qmax)));
[~,kDGR] = max(DGR);
%kGR=kGR-1;

% compute  DDR and kDDR
den = max((mu(2:(qmax+1))-mu(3:(qmax+2))),mum);
ncorrections = sum(den == mum);
DDR = (mu(1:(qmax))-mu(2:(qmax+1)))./den;
[~,kDDR] = max(DDR);

