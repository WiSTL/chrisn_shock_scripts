function [f,A,c] = C_pdf(C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

I = round(C,3);

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

A = interp1(c,A,[0.001:0.001:1-0.001]);
c=[0.001:0.001:1-0.001];

dc = diff(smooth(c,25));
dA = diff(smooth(A,25))';

[nx,ny] = size(I);

f = smooth(abs((dA./dc')/(nx*ny)),25);%Probability Density Function

end

