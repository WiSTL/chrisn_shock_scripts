function w = chris_fourier_derivative(f,a,b,n)
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
Nx = length(f);

% fftx = fft_corrected(f',1);
fftx = fft(f);
k = 2*pi/((b-a))*[0:Nx/2-1, 0, -Nx/2+1:-1]';

if length(fftx)~=length(k)
    dffft = ((1i*k).^n).*fftx(2:end);
else
    dffft = ((1i*k).^n).*fftx;
end
w = real(ifft(dffft));

end