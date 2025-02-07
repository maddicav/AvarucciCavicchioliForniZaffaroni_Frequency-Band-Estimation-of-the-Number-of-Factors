%--------------------------------------------------------------------------
% This file load the data from "FREDMD.xlxs" and transfom them 
% according to the transformation codes
% 1 = no transformation
% 2 = first difference
% 4 = log transformation
% 5 = log difference
% 7 = first difference of the percentage variation
%--------------------------------------------------------------------------
ldata=xlsread('FRED-MD.xlsx','current_MD','B15:DS739'); 
trans=xlsread('FRED-MD.xlsx','current_MD','B4:DS4');
for i=1:length(trans)
    if trans(i)==1
        x(:,i)=ldata(3:end,i);
    elseif trans(i)==2
        x(:,i)=diff(ldata(2:end,i));
    elseif trans(i)==4
        x(:,i)=log(ldata(3:end,i))*100;
    elseif trans(i)==5
        x(:,i)=diff(log(ldata(2:end,i)))*100;
    elseif trans(i)==7
        x(:,i)=diff(ldata(2:end,i)./ldata(1:end-1,i)-1);
    end
end

