function [dfdx, dfdy] = chris_fourier_gradient(f,dx)
% FOURIERDERIVATIVE Fourier derivative
%     dfdx = FOURIERDERIVATIVE(f,a,b) approximates the derivative a
%     discrete function f over the domain (a,b).  f, a vector, must be
%     uniformly sampled, periodic, and contain an even number of samples.
%     For best results, f should be periodic such that f(x + a) = f(x + b).
%     As an example,
%
%          x = linspace(0,pi);
%          f = exp(cos(x).*sin(2*x));
%          dfdx = fourierderivative(f,0,pi);
%
%     Results for nonperiodic f are dubious.

[Ny,Nx] = size(f);
[kx, ky] = meshgrid((-Nx/2:Nx/2)/Nx*dx, (-Ny/2:Ny/2)/Ny*dx);

[nky,nkx] = size(kx)
[nfy,nfx] = size(fft2(f))

if nky == nfy
    dfdx = real(ifft2(1i .* kx .* fft2(f)));
    dfdy = real(ifft2(1i .* ky .* fft2(f)));
else
    dfdx = real(ifft2(1i .* kx(2:end,2:end) .* fft2(f)));
    dfdy = real(ifft2(1i .* ky(2:end,2:end) .* fft2(f)));
end



% Nx = length(f);
% 
% % fftx = fft_corrected(f',1);
% fftx = fft(f);
% k = 2*pi/((b-a))*[0:Nx/2-1, 0, -Nx/2+1:-1]';
% 
% if length(fftx)~=length(k)
%     dffft = ((1i*k).^n).*fftx(2:end);
% else
%     dffft = ((1i*k).^n).*fftx;
% end
% w = real(ifft(dffft));

end