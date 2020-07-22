function [out,out2] = Alex_notch_filter(I,A,B)

%%%%
%
% A = 20, B = 40
%
%%%%

[origin, lines] = laserorigin(I);

out2 = I;%lasertransform(I, origin, 'forward');

[ni, nj] = size(out2);
% Select regions to remove: horizontal stripes corresponding to bands
rm_i = (floor(ni/2) - A):(ceil(ni/2) + A);
rm_j = [(floor(nj/4):floor(nj/2))-B (floor(nj/2):3*floor(nj/4))+B];

F = fftshift(fft2(out2)); % FFT to freq. domain, swap quadrants for clarity
F(rm_i, rm_j) = 0;       % Remove banded regions
M = (gausswin(ni,1) * gausswin(nj, 1)'); % Use Gaussian to remove noise
G = M .* F;
out = real(ifft2(ifftshift(G)));

out = out;%lasertransform(out, origin, 'inverse');

end

