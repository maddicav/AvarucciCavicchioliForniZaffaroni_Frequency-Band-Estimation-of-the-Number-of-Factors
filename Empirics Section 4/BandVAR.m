
function [idirfs,idbootirfs, idshocks, D] = BandVAR(irfs, bootirfs, shocks, bnd, kk, T, idsign, kkk,SS)
%% variable number of inputs
if nargin < 9
    
    ss = sqrt(sum(sum(irfs.^2,3),2));
else
    
ss= sqrt(SS);
end
if nargin < 8
    kkk=1;
end
if nargin < 7
    idsign = ones(1,size(shocks,2))*3;
end


if nargin < 6
   T = 256;
end
%% band of interest

%int = pi/T;
int = 2*pi/T;
omega = int:int:pi;
ome = omega(omega > bnd(1) & omega < bnd(2));



%%  point estimates

[idirfs, idshocks, D] = pointestimate(irfs, shocks,ome,kk,idsign,kkk,ss);

 %%  boot 
 
 if isempty(bootirfs) 
     idbootirfs = [];
 else
     idbootirfs = bootirfs*nan;
   for j = 1:size(bootirfs,4)
      idbootirfs(:,:,:,j) = pointestimate(bootirfs(:,:,:,j), shocks,ome,kk,idsign,kkk,ss); 
   end
 end
 

function [idirfs, idshocks, D] = pointestimate(irfs, shocks, ome, kk, idsign, kkk,ss)

ns = size(irfs,2);
nom = length(ome);
h = size(irfs,3);





ck = nan(length(kk), ns,nom); ckck=nan(ns,ns,nom);
for j=1:nom
    for jj=1:length(kk)
        ck(jj,:,j) = (1/ss(kk(jj)))*squeeze(irfs(kk(jj),:,:))*exp(-1i*ome(j)*(0:h-1)');
    end
    ckck(:,:,j) = ck(:,:,j)'*ck(:,:,j);
end




thk = mean(real(ckck),3);

% compute the egenvectors P
[P,D] = eig(thk); 
[D,o] = sort(diag(D),'descend' ); 
P = P(:,o);
%check sign
for i = 1:length(idsign)
P(:,i) = P(:,i)*sign(sum(irfs(kkk,:,1:idsign(i)),3)*P(:,i));
end
% compute the identified irfs and shocks
idirfs  = polynomialmatricesproduct(irfs,P,h);
idshocks = shocks*P;