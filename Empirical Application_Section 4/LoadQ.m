
ldata=xlsread('FREDreduced.xlsx','Foglio1','B6:HI249'); 
trans=xlsread('FREDreduced.xlsx','Foglio1','B4:HI4');
trans1=xlsread('FREDreduced.xlsx','Foglio1','B5:HI5');
%trans=trans1;ldata=ldata(2:end,:);
for i=1:length(trans)
    if trans(i)==1
        x(:,i)=ldata(3:end,i);
    elseif trans(i)==2
        x(:,i)=diff(ldata(2:end,i));
     elseif trans(i)==3
         x(:,i)=diff(diff(ldata(1:end,i)));
    elseif trans(i)==4
        x(:,i)=log(ldata(3:end,i))*100;
    elseif trans(i)==5
        x(:,i)=diff(log(ldata(2:end,i)))*100;
     elseif trans(i)==6
        %x(:,i)=diff(log(ldata(2:end,i)))*100;
        x(:,i)=diff(diff(log(ldata(1:end,i))))*100;
    elseif trans(i)==7
        x(:,i)=diff(ldata(2:end,i)./ldata(1:end-1,i)-1);
    end
end

