function [B] = chris_pcolor(A,range,y)

B = squeeze(A);

if y==1
figure
end
pcolor(B),colorbar,axis equal, axis off, shading flat;
caxis(range);


end

