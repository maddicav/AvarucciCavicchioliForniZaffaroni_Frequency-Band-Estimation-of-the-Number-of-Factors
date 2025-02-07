sy=245; % number of other static endogenous variables
sx=90; % number of predetermined variables
sp=91;  % number of jump variables
sz=33;  % number of exogenous shocks
J=30;   % number of sectors

% Next, define some auxiliary matrices

Omega21=Omega2(:,1:sx);
Omega22=Omega2(:,sx+1:end);
M1=A11+A12*F1;
M2=B1+A12*F2;
M3=inv(eye(sy)-Omega1)*Omega21;
M4=inv(eye(sy)-Omega1)*Omega22;
M5=inv(eye(sy)-Omega1)*Omega3;

% Finally, define the matrices selecting the variables in our dataset from the 
% set of endogenous static variables, the set of the predetermined 
% variables and the set of the jump variables

S1=[eye(2*J) zeros(2*J,sy-2*J);...
        zeros(J,3*J) eye(J) zeros(J,sy-4*J);...
        zeros(J,7*J) eye(J) zeros(J,sy-8*J);...
        zeros(4,8*J) eye(4) zeros(4,sy-4-8*J)];
S2=[zeros(1,J) 1 zeros(1,sx-J-1)];
S3=[zeros(J,J-1) 4*eye(J) zeros(J,sp+1-2*J);...
        zeros(1,2*J-1) 4 zeros(1,sp-2*J)];