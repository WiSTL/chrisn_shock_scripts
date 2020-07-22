function [lambdacx,lambdacy] = C_length_scales(A)

dx = 0.0315e-2;
dy = dx;

% [~,A,~] = C_fluctuation_595(A);

[nx,ny] = size(A);

xx=1:nx;
yy=1:ny;

M = length(xx);
N = length(yy);
k_x = (0:M-1)/dx/M;
k_y = (0:N-1)/dx/N;

cmx = mean(A');

difcm = mean(diff(smooth(cmx,25)')./dy).^2;
cmm=mean(cmx).^2;

lambdacx = (sqrt(cmm/difcm));

cmy = mean(A);
cmm=mean(cmy).^2;

difcm = mean(diff(smooth(cmy,25)')./dx).^2;

lambdacy = (sqrt(cmm/difcm));


end

