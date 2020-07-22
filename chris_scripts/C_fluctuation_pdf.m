function [Ar,c] = C_fluctuation_pdf(C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,Cm, meanC] = C_fluctuation_595(C);

I = round(sqrt(Cm.^2),3);

A=[];
n=0;

for c=0.001:0.001:1-0.001
    n=n+1;
    [A(n),~] = size(I(I==c));    
end

c=[0.001:0.001:1-0.001];


for i=1:n
   if(A(i)==0)
       A(i) = nan;
       c(i) = nan;
   end
end

A=A(~any(isnan(A),1));%Area for a given Concentration
c=c(~any(isnan(c),1));

Ar = interp1(c,A,[0:0.001:0.5]);
c=[0:0.001:0.5];


end

