%--------------------------------------------------------------------------
% This file load the data from "FREDreduced.xlxs" and transfom them 
% according to the transformation codes
% 1 = no transformation
% 2 = first difference
% 5 = log diffeence
% 7 = first difference of the percentage variation
%--------------------------------------------------------------------------
ldata=xlsread('FREDreduced.xlsx','Foglio1','B5:HI248'); % 216 variables
trans=xlsread('FREDreduced.xlsx','Foglio1','B4:HI4');
for i=1:length(trans)
    if trans(i)==1
        x(:,i)=ldata(3:end,i);
    elseif trans(i)==2
        x(:,i)=diff(ldata(2:end,i));
    elseif trans(i)==5
        x(:,i)=diff(log(ldata(2:end,i)))*100;
    elseif trans(i)==7
        x(:,i)=diff(ldata(2:end,i)./ldata(1:end-1,i)-1);
    end
end

