function [meanI] = shear_layer_mean(I)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[nx ny ni] = size(I)

meanI=zeros(551,ny);

for i = 1:ni
   [~,~,~, i5, i95] = C_shannon_entropy_595(I(:,:,i));
   i50 = 0.5*(i5+i95);
   meanI = (meanI+I(i50-250:i50+300,:,i)); 
end

meanI = meanI/ni;

end

