function [s_xt,s_yt,X,Y,T] = chris_phantom_xt(fn)

[ I,~] = readcine(fn);
[nx,ny,ni] = size(I);
Ibg = mean(I(:,:,end-20:end),3);

dt = 5e-5;
ni = 500;

xlin=linspace(4.6645,4.7661,ny);
ylin=linspace(-0.3175,0.3175,nx);
tlin=linspace(0.01,dt*ni+0.01,ni);

for i = 1:ni
Is(:,:,i) = I(:,:,i) - Ibg;
end

for i = 1:ni
s_yt(:,i) = mean(Is(:,:,i)');
end
for i = 1:ni
s_xt(:,i) = mean(Is(:,:,i));
end

[X,T]=meshgrid(xlin,tlin);
figure
pcolor(X,T,flipud(s_xt)'), shading flat
colormap(flipud(bone));
caxis([0 1.1256e+03])
title(fn)

s_xt = flipud(s_xt)';

[Y,~]=meshgrid(ylin,tlin);
% figure(2)
% pcolor(Y,T,flipud(s_yt)'), shading flat
% colormap(flipud(bone));


end

