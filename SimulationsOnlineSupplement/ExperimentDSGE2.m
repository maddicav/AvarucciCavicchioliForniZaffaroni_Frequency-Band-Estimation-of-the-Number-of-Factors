%-----------------------------------------------------------------------------------
% Output: TABLE S.8 in Avarucci et al 2024, Online Supplement 
%
% Top Panel (n=156, T=120): select data_generator in Line 16  
% Bottom Panel (n=156, T=120): select data_generator1 in Line 17
%
% The files data_generator  and data_generator1 generate the data according
% to A. Onatski and F. Ruge-Murcia (2013)"Factor Analysis of a Large DSGE
% Model" and are available here  https://sites.google.com/site/frugemurcia/welcome/replication-files
%------------------------------------------------------------------------------------
clear all
nreps = 500;
%
for j = 1:nreps
    j
% data_generator
data_generator1;
x = standardize(DATA');
[kDER(j), kDGR(j), kDDR(j)] = ACFZcrit(x, 17,.75);
end
a = [kDER' kDGR' kDDR' ];
a1 = sum(a==1);
a2 = sum(a==2);
a3 = sum(a==3);
a4 = sum(a>=4);
b=[a1;a2;a3;a4]'/nreps*100;
display('Table S.8')
disp(b)