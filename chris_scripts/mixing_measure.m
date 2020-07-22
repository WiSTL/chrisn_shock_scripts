function [d] = mixing_measure(C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,ny] = size(C);

d = [];

for i = 1:ny
    cr = C(:,i);
    
    l = find(cr == 0);
    
    l1 = [l;0];
    l2 =[0;l];
    
    d(i) = max(abs(l1-l2));
end

d=min(d)/nanmax(d);

end

