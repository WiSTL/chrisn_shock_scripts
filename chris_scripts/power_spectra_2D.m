function [S2d, S2dfilt, k_x, k_y, S2dmirror, kxfull, kyfull] = power_spectra_2D(C,dx,dy)

[nx,ny] = size(C);

xx=1:nx;
yy=1:ny;

M = length(xx);
N = length(yy);
k_x = (0:M-1)/dx/M;
k_y = (0:N-1)/dy/N;

wf2d=(1/M/dx)*(1/N/dy)*fft2(C);
S2d=(wf2d.*conj(wf2d));

have=fspecial('average',[7 7]);
S2dfilt=imfilter(S2d,have,'replicate');

t=k_x(2:M/2);
t1=-fliplr(t);
kxfull=[t1,t];
t=k_y(2:N/2);
t1=-fliplr(t);
kyfull=[t1,t];
est=S2dfilt(2:M/2,2:N/2);
est1=fliplr(est);
est2=[est1,est];
est3=flipud(est2);
S2dmirror=[est3;est2];

end

