function [A] = mixing_region_area(I)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[nx, ny] = size(I);

I_1 = I;
I_10 = I;


countzero = length(I(I<=0.05));
countone = length(I(I>=0.95));


for i=1:nx
    for j=1:ny
        I_1(i,j) = (I_1(i,j)>= 0.49);
    end
end

for i=1:nx
    for j=1:ny
        I_10(i,j) = (I_10(i,j)>= 0.51);
    end
end

% for i=1:nx
%     for j=1:ny
%         I_50(i,j) = (I_50(i,j)== 0.5);
%     end
% end
% 
% for i=1:nx
%     for j=1:ny
%         I_45(i,j) = (I_45(i,j)== 0.45);
%     end
% end
% 
% for i=1:nx
%     for j=1:ny
%         I_55(i,j) = (I_55(i,j)== 0.55);
%     end
% end

% A = (nx*ny - (countzero+countone))/nx;
A = sum(sum(I_1-I_10))/nx;

end

