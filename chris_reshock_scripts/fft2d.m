function [b,kx,ky] = fft2d(a)

[ny,nx] = size(a);

kx=(2*pi)*(1:nx)/(nx);
ky=(2*pi)*(1:ny)/(ny);

b = fft(a, [], 1);
b = fft(b, [], 2);
end

