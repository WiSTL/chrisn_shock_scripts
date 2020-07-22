function IMAGEStoGIF(A,range,filename,dt,cmap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nx, ny, nimg] = size(A);

h = figure('Position', [100, 100, 800, 800]);

pcolor(flipud(squeeze(A(:,:,1)))), shading flat
caxis(range)
set(gca,'Ydir','reverse');
set(gcf,'color','w');
% set(gca,'xscale','log')
axis equal
axis tight manual 
axis off
% colormap gray
cmocean(cmap)
% colorbar
text(20,5,strcat(['t: ',num2str(1)]),'FontSize',18)
set(gca,'nextplot','replacechildren');

f=getframe(h);
[im,cm] = rgb2ind(f.cdata,256);

imwrite(im,cm,filename,'DelayTime',dt,'LoopCount',inf);

for i = 2:nimg

pcolor(flipud(squeeze(A(:,:,i)))), shading flat
caxis(range)
set(gca,'Ydir','reverse');
set(gcf,'color','w');
% set(gca,'xscale','log')
axis equal
axis off
% colormap gray
cmocean(cmap)
% colorbar
text(20,5,strcat(['t: ',num2str(i)]),'FontSize',18)

f=getframe(h);
[im,cm] = rgb2ind(f.cdata,256);

imwrite(im,cm,filename,'DelayTime',dt,'writemode','append');

end

close(h);

end
