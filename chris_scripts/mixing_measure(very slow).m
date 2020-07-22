function [m] = mixing_measure(C)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

[nx,ny] = size(C);

C2 = round(1-C,4);

d=[];

n=0;

for i = 1:nx
    for j = 1:ny
        c = round(C(i,j),4);
        [row,col]=find(C2==c);
        
        if(isempty(row) ~= 1)
            n=n+1;
            a=sqrt((i-row).^2 + (j-col).^2);
            d(n) = min(a);
        end
    end
end

m = sum(d.^2)/(length(d));

end

