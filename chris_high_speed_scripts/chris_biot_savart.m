function [u,v] = chris_biot_savart(omega,dx,dy)

[~,nxx] = size(omega);

omega = [fliplr(omega) omega fliplr(omega)];

[ny,nx] = size(omega);

x = linspace(-dx*nxx/2,dx*nxx/2,nxx);
y = linspace(-dy*ny/2,dy*ny/2,ny);

[X,Y] = meshgrid(x,y);

R = sqrt(X.^2+Y.^2);

Kx = -Y./R.^2;
Ky = X./R.^2;

Kx(isnan(Kx)) = 0;
Ky(isnan(Ky)) = 0;

u = imfilter(omega,Kx/(2*pi),'conv');
v = imfilter(omega,Ky/(2*pi),'conv');

u = u(:,nxx+1:2*nxx);
v = v(:,nxx+1:2*nxx);

end

