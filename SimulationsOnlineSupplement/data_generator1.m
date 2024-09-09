% This file generates data using a backward looking dynamic system
% which is the solution to the forward looking equations of BCR model
loadings
parameters

% innovations are standard normals
IN=randn(33,242);
% initialize by a draw from ergodic distribution
SIG=diag(sigma.*(ones(33,1)-rho.^2).^(-1/2));
aZ=SIG*IN(:,1);
VX=dlyap(M1,M2*SIG*SIG'*M2');
SVX=chol((VX+VX')/2+0.00001*eye(sx));
aX=SVX'*randn(sx,1);
RHO=diag(rho);
% the data simulation cycle 
for t=1:241
    aX(:,t+1)=M1*aX(:,t)+M2*aZ(:,t);
    aP(:,t)=F1*aX(:,t)+F2*aZ(:,t);
    aY(:,t)=M3*aX(:,t)+M4*aP(:,t)+M5*aZ(:,t);
    aZ(:,t+1)=RHO*aZ(:,t)+diag(sigma)*IN(:,t+1);
end
% selection of relevant variables
DATA=[S1*aY;S2*aX(:,2:end);S3*aP];
DATA=[DATA(1:end-1-J,2:end);...
        DATA(end-J:end-1,2:end)-DATA(end-J:end-1,1:end-1)+kron(ones(J,1),DATA(end,2:end));...
        DATA(end,2:end)];