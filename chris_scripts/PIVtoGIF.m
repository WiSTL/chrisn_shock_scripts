function PIVtoGIF(A,B,range,filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nx, ny] = size(A);

h = figure('Position', [100, 100, ny, nx]);

imagesc(A), shading flat
caxis([0 150])
axis equal
axis off
colormap gray;

t = text(50,50,'frame 1','color','w','FontSize',20);


f=getframe;
[im,cm] = rgb2ind(f.cdata,256);

imwrite(im,cm,filename,'DelayTime',0.5,'LoopCount',inf);

imagesc(B), shading flat
caxis([0 200])
axis equal
axis off

t = text(50,50,'frame 2','color','w','FontSize',20);


f=getframe;
[im, cm] = rgb2ind(f.cdata,256);

imwrite(im,cm,filename,'DelayTime',0.5,'writemode','append');

close(h);

end

